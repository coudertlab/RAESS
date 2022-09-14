// Copyright 2020 Global Phasing Ltd.
//
// Class Intensities that reads multi-record data from MTZ, mmCIF or XDS_ASCII
// and merges it into mean or anomalous intensities.
// It can also read merged data.

#ifndef GEMMI_MERGE_HPP_
#define GEMMI_MERGE_HPP_

#include <cassert>
#include "symmetry.hpp"
#include "unitcell.hpp"
#include "util.hpp"     // for vector_remove_if
#include "mtz.hpp"      // for Mtz
#include "refln.hpp"    // for ReflnBlock
#include "xds_ascii.hpp" // for XdsAscii

namespace gemmi {

inline std::string miller_str(const Miller& hkl) {
  std::string s;
  for (int i = 0; i < 3; ++i) {
    s += (i == 0 ? '(' : ' ');
    s += std::to_string(hkl[i]);
  }
  s += ')';
  return s;
}

struct Intensities {
  enum class Type { None, Unmerged, Mean, Anomalous };

  struct Refl {
    Miller hkl;
    int isign;  // 1 for I(+), -1 for I(-), 0 for mean
    double value;
    double sigma;

    bool operator<(const Refl& o) const {
      return std::tie(hkl[0], hkl[1], hkl[2], isign) <
             std::tie(o.hkl[0], o.hkl[1], o.hkl[2], o.isign);
    }
  };

  std::vector<Refl> data;
  const SpaceGroup* spacegroup = nullptr;
  UnitCell unit_cell;
  double wavelength;
  Type type = Type::None;

  static const char* type_str(Type itype) {
    switch (itype) {
      case Type::None: return "n/a";
      case Type::Unmerged: return "I";
      case Type::Mean: return "<I>";
      case Type::Anomalous: return "I+/I-";
    }
    unreachable();
  }

  const char* type_str() const { return Intensities::type_str(type); }

  std::string spacegroup_str() const { return spacegroup ? spacegroup->xhm() : "none"; }

  std::array<double,2> resolution_range() const {
    double min_1_d2 = INFINITY;
    double max_1_d2 = 0;
    for (const Refl& x : data) {
      double a_1_d2 = unit_cell.calculate_1_d2(x.hkl);
      if (a_1_d2 < min_1_d2)
        min_1_d2 = a_1_d2;
      if (a_1_d2 > max_1_d2)
        max_1_d2 = a_1_d2;
    }
    return {{ 1 / std::sqrt(min_1_d2), 1 / std::sqrt(max_1_d2) }};
  }

  // pre: both are sorted
  Correlation calculate_correlation(const Intensities& other) const {
    Correlation corr;
    auto r1 = data.begin();
    auto r2 = other.data.begin();
    while (r1 != data.end() && r2 != other.data.end()) {
      if (r1->hkl == r2->hkl && r1->isign == r2->isign) {
        corr.add_point(r1->value, r2->value);
        ++r1;
        ++r2;
      } else if (*r1 < *r2) {
        ++r1;
      } else {
        ++r2;
      }
    }
    return corr;
  }

  void remove_systematic_absences() {
    if (!spacegroup)
      return;
    GroupOps gops = spacegroup->operations();
    vector_remove_if(data, [&](Refl& x) {
        return gops.is_systematically_absent(x.hkl);
    });
  }

  void sort() { std::sort(data.begin(), data.end()); }

  void merge_in_place(Type itype) {
    if (data.empty())
      return;
    if (itype == Type::Mean)
      // discard signs so that merging produces Imean
      for (Refl& refl : data)
        refl.isign = 0;
    sort();
    std::vector<Refl>::iterator out = data.begin();
    double sum_wI = 0.;
    double sum_w = 0.;
    for (auto in = data.begin(); in != data.end(); ++in) {
      if (out->hkl != in->hkl || out->isign != in->isign) {
        out->value = sum_wI / sum_w;
        out->sigma = 1.0 / std::sqrt(sum_w);
        sum_wI = sum_w = 0.;
        ++out;
        out->hkl = in->hkl;
        out->isign = in->isign;
      }
      double w = 1. / (in->sigma * in->sigma);
      sum_wI += w * in->value;
      sum_w += w;
    }
    out->value = sum_wI / sum_w;
    out->sigma = 1.0 / std::sqrt(sum_w);
    data.erase(++out, data.end());
  }

  // for unmerged centric reflections set isign=1.
  void switch_to_asu_indices(bool merged=false) {
    GroupOps gops = spacegroup->operations();
    ReciprocalAsu asu(spacegroup);
    for (Refl& refl : data) {
      if (asu.is_in(refl.hkl)) {
        if (!merged && refl.isign == -1 && gops.is_reflection_centric(refl.hkl))
          refl.isign = 1;
        continue;
      }
      auto hkl_isym = asu.to_asu(refl.hkl, gops);
      refl.hkl = hkl_isym.first;
      if (!merged) {
        if (gops.is_reflection_centric(refl.hkl))
          refl.isign = 1;
        else
          refl.isign = (hkl_isym.second % 2 == 0 ? -1 : 1);
      }
    }
  }

  void read_unmerged_intensities_from_mtz(const Mtz& mtz) {
    if (mtz.batches.empty())
      fail("expected unmerged file");
    const Mtz::Column* isym_col = mtz.column_with_label("M/ISYM");
    if (!isym_col || isym_col->idx != 3)
      fail("unmerged file should have M/ISYM as 4th column");
    const Mtz::Column& col = mtz.get_column_with_label("I");
    size_t value_idx = col.idx;
    size_t sigma_idx = mtz.get_column_with_label("SIGI").idx;
    unit_cell = mtz.get_average_cell_from_batch_headers(nullptr);
    spacegroup = mtz.spacegroup;
    if (!spacegroup)
      fail("unknown space group");
    wavelength = mtz.dataset(col.dataset_id).wavelength;
    for (size_t i = 0; i < mtz.data.size(); i += mtz.columns.size()) {
      int isign = ((int)mtz.data[i + 3] % 2 == 0 ? -1 : 1);
      add_if_valid(mtz.get_hkl(i), isign, mtz.data[i + value_idx], mtz.data[i + sigma_idx]);
    }
    // Aimless >=0.7.6 (from 2021) has an option to output unmerged file
    // with original indices instead of reduced indices, with all ISYM = 1.
    switch_to_asu_indices();
    type = Type::Unmerged;
  }

  void read_mean_intensities_from_mtz(const Mtz& mtz) {
    if (!mtz.batches.empty())
      fail("expected merged file");
    const Mtz::Column* col = mtz.column_with_one_of_labels({"IMEAN", "I"});
    if (!col)
      fail("Mean intensities (IMEAN or I) not found.");
    size_t sigma_idx = mtz.get_column_with_label("SIG" + col->label).idx;
    copy_metadata(mtz);
    wavelength = mtz.dataset(col->dataset_id).wavelength;
    read_data(MtzDataProxy{mtz}, col->idx, sigma_idx);
    type = Type::Mean;
  }

  void read_anomalous_intensities_from_mtz(const Mtz& mtz, bool check_complete=false) {
    if (!mtz.batches.empty())
      fail("expected merged file");
    const Mtz::Column& col = mtz.get_column_with_label("I(+)");
    size_t value_idx[2] = {col.idx, mtz.get_column_with_label("I(-)").idx};
    size_t sigma_idx[2] = {mtz.get_column_with_label("SIGI(+)").idx,
                           mtz.get_column_with_label("SIGI(-)").idx};
    int mean_idx = -1;
    if (check_complete)
      if (const Mtz::Column* mean_col = mtz.column_with_one_of_labels({"IMEAN", "I"}))
        mean_idx = (int) mean_col->idx;
    copy_metadata(mtz);
    wavelength = mtz.dataset(col.dataset_id).wavelength;
    read_anomalous_data(MtzDataProxy{mtz}, mean_idx, value_idx, sigma_idx);
    type = Type::Anomalous;
  }

  void read_merged_intensities_from_mtz(const Mtz& mtz) {
    if (mtz.column_with_label("I(+)"))
      read_anomalous_intensities_from_mtz(mtz, true);
    else
      read_mean_intensities_from_mtz(mtz);
  }

  void read_mtz(const Mtz& mtz, Type itype) {
    switch (itype) {
      case Type::Unmerged:
        read_unmerged_intensities_from_mtz(mtz);
        break;
      case Type::Mean:
        read_mean_intensities_from_mtz(mtz);
        break;
      case Type::Anomalous:
        read_anomalous_intensities_from_mtz(mtz);
        break;
      case Type::None:
        assert(0);
        break;
    }
  }

  void read_unmerged_intensities_from_mmcif(const ReflnBlock& rb) {
    size_t value_idx = rb.get_column_index("intensity_net");
    size_t sigma_idx = rb.get_column_index("intensity_sigma");
    copy_metadata(rb);
    wavelength = rb.wavelength;
    read_data(ReflnDataProxy(rb), value_idx, sigma_idx);
    switch_to_asu_indices();
    type = Type::Unmerged;
  }

  void read_mean_intensities_from_mmcif(const ReflnBlock& rb) {
    size_t value_idx = rb.get_column_index("intensity_meas");
    size_t sigma_idx = rb.get_column_index("intensity_sigma");
    copy_metadata(rb);
    wavelength = rb.wavelength;
    read_data(ReflnDataProxy(rb), value_idx, sigma_idx);
    type = Type::Mean;
  }

  void read_anomalous_intensities_from_mmcif(const ReflnBlock& rb,
                                             bool check_complete=false) {
    size_t value_idx[2] = {rb.get_column_index("pdbx_I_plus"),
                           rb.get_column_index("pdbx_I_minus")};
    size_t sigma_idx[2] = {rb.get_column_index("pdbx_I_plus_sigma"),
                           rb.get_column_index("pdbx_I_minus_sigma")};
    int mean_idx = -1;
    if (check_complete)
      mean_idx = rb.find_column_index("intensity_meas");
    copy_metadata(rb);
    wavelength = rb.wavelength;
    read_anomalous_data(ReflnDataProxy(rb), mean_idx, value_idx, sigma_idx);
    type = Type::Anomalous;
  }

  void read_merged_intensities_from_mmcif(const ReflnBlock& rb) {
    if (rb.find_column_index("pdbx_I_plus"))
      read_anomalous_intensities_from_mmcif(rb, true);
    else
      read_mean_intensities_from_mmcif(rb);
  }

  void read_mmcif(const ReflnBlock& rb, Type itype) {
    switch (itype) {
      case Type::Unmerged:
        read_unmerged_intensities_from_mmcif(rb);
        break;
      case Type::Mean:
        read_mean_intensities_from_mmcif(rb);
        break;
      case Type::Anomalous:
        read_anomalous_intensities_from_mmcif(rb);
        break;
      case Type::None:
        break;
    }
  }

  void read_unmerged_intensities_from_xds(const XdsAscii& xds) {
    unit_cell = xds.unit_cell;
    spacegroup = find_spacegroup_by_number(xds.spacegroup_number);
    wavelength = xds.wavelength;
    data.reserve(xds.data.size());
    for (const XdsAscii::Refl& in : xds.data)
      add_if_valid(in.hkl, 0, in.iobs, in.sigma);
    switch_to_asu_indices();
    type = Type::Unmerged;
  }

private:
  template<typename Source>
  void copy_metadata(const Source& source) {
    unit_cell = source.cell;
    spacegroup = source.spacegroup;
    if (!spacegroup)
      fail("unknown space group");
  }

  void add_if_valid(const Miller& hkl, int isign, double value, double sigma) {
    // XDS marks rejected reflections with negative sigma.
    // Sigma 0.0 is also problematic - it rarely happens (e.g. 5tkn).
    if (!std::isnan(value) && sigma > 0)
      data.push_back({hkl, isign, value, sigma});
  }

  template<typename DataProxy>
  void read_data(const DataProxy& proxy, size_t value_idx, size_t sigma_idx) {
    for (size_t i = 0; i < proxy.size(); i += proxy.stride())
      add_if_valid(proxy.get_hkl(i), 0,
                   proxy.get_num(i + value_idx),
                   proxy.get_num(i + sigma_idx));
  }

  template<typename DataProxy>
  void read_anomalous_data(const DataProxy& proxy, int mean_idx,
                           size_t (&value_idx)[2], size_t (&sigma_idx)[2]) {
    GroupOps gops = spacegroup->operations();
    for (size_t i = 0; i < proxy.size(); i += proxy.stride()) {
      Miller hkl = proxy.get_hkl(i);
      bool centric = gops.is_reflection_centric(hkl);
      if (mean_idx >= 0 && !std::isnan(proxy.get_num(i + mean_idx)) && !centric) {
        if (std::isnan(proxy.get_num(i + value_idx[0])) &&
            std::isnan(proxy.get_num(i + value_idx[1])))
          fail(miller_str(hkl), " has <I>, but I(+) and I(-) are both null");
      }
      add_if_valid(hkl, 1, proxy.get_num(i + value_idx[0]),
                           proxy.get_num(i + sigma_idx[0]));
      if (!centric)  // ignore I(-) of centric reflections
        add_if_valid(hkl, -1, proxy.get_num(i + value_idx[1]),
                              proxy.get_num(i + sigma_idx[1]));
    }
  }
};

} // namespace gemmi
#endif
