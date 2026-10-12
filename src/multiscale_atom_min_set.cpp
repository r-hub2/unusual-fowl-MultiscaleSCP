// MultiscaleSCP-owned model storage. Data crosses the package boundary as an
// R list; no prioritizr C++ object is cast or dereferenced here.

#include <Rcpp.h>
using namespace Rcpp;

#include <vector>
#include <string>
#include <algorithm>

class MultiscaleModel {
public:
  std::string _modelsense;
  std::size_t _number_of_features;
  std::size_t _number_of_planning_units;
  std::size_t _number_of_zones;
  std::vector<std::size_t> _A_i;
  std::vector<std::size_t> _A_j;
  std::vector<double> _A_x;
  std::vector<double> _obj;
  std::vector<double> _lb;
  std::vector<double> _ub;
  std::vector<double> _rhs;
  std::vector<std::string> _sense;
  std::vector<std::string> _vtype;
  std::vector<std::string> _row_ids;
  std::vector<std::string> _col_ids;
  bool _compressed_formulation;
};

static std::size_t max_row_index_plus_one(const std::vector<std::size_t>& rows) {
  if (rows.empty()) return 0;
  return (*std::max_element(rows.begin(), rows.end())) + 1;
}

static void ensure_col_vectors(MultiscaleModel* ptr, std::size_t ncol) {
  if (ptr->_obj.size() < ncol) ptr->_obj.resize(ncol, 0.0);
  if (ptr->_lb.size() < ncol) ptr->_lb.resize(ncol, 0.0);
  if (ptr->_ub.size() < ncol) ptr->_ub.resize(ncol, 1.0);
  if (ptr->_vtype.size() < ncol) ptr->_vtype.resize(ncol, "B");
  if (ptr->_col_ids.size() < ncol) ptr->_col_ids.resize(ncol, "pu");
}

static void add_neutral_rows_until(MultiscaleModel* ptr, std::size_t nrow) {
  while (ptr->_rhs.size() < nrow) {
    ptr->_rhs.push_back(0.0);
    ptr->_sense.push_back(">=");
    ptr->_row_ids.push_back("neutral_rij_target");
  }
}

// [[Rcpp::export]]
Rcpp::List rcpp_build_multiscale_atom_min_set_model(
  const Rcpp::List model,
  const Rcpp::List targets_list,
  const Rcpp::NumericVector atom_costs,
  const Rcpp::List cover_by_atom,
  SEXP feature_atom_matrixSEXP
) {
  MultiscaleModel storage;
  MultiscaleModel* ptr = &storage;
  ptr->_modelsense = Rcpp::as<std::string>(model["modelsense"]);
  ptr->_number_of_features = Rcpp::as<std::size_t>(model["number_of_features"]);
  ptr->_number_of_planning_units = Rcpp::as<std::size_t>(model["number_of_planning_units"]);
  ptr->_number_of_zones = Rcpp::as<std::size_t>(model["number_of_zones"]);
  ptr->_A_i = Rcpp::as<std::vector<std::size_t>>(model["A_i"]);
  ptr->_A_j = Rcpp::as<std::vector<std::size_t>>(model["A_j"]);
  ptr->_A_x = Rcpp::as<std::vector<double>>(model["A_x"]);
  ptr->_obj = Rcpp::as<std::vector<double>>(model["obj"]);
  ptr->_lb = Rcpp::as<std::vector<double>>(model["lb"]);
  ptr->_ub = Rcpp::as<std::vector<double>>(model["ub"]);
  ptr->_rhs = Rcpp::as<std::vector<double>>(model["rhs"]);
  ptr->_sense = Rcpp::as<std::vector<std::string>>(model["sense"]);
  ptr->_vtype = Rcpp::as<std::vector<std::string>>(model["vtype"]);
  ptr->_row_ids = Rcpp::as<std::vector<std::string>>(model["row_ids"]);
  ptr->_col_ids = Rcpp::as<std::vector<std::string>>(model["col_ids"]);
  ptr->_compressed_formulation = Rcpp::as<bool>(model["compressed_formulation"]);


  if (ptr->_number_of_zones != 1) {
    Rcpp::stop("Multiscale atom min-set currently supports single-zone problems only.");
  }

  if (!ptr->_compressed_formulation) {
    Rcpp::stop("Multiscale atom min-set currently supports the compressed formulation only.");
  }

  // Read sparse Matrix::dgCMatrix directly without Armadillo.
  Rcpp::S4 feature_atom_matrix(feature_atom_matrixSEXP);

  Rcpp::IntegerVector dims = feature_atom_matrix.slot("Dim");
  Rcpp::IntegerVector i_slot = feature_atom_matrix.slot("i");
  Rcpp::IntegerVector p_slot = feature_atom_matrix.slot("p");
  Rcpp::NumericVector x_slot = feature_atom_matrix.slot("x");

  if (dims.size() != 2) {
    Rcpp::stop("feature_atom_matrix must be a two-dimensional sparse Matrix object.");
  }

  const std::size_t n_pu =
    ptr->_number_of_planning_units;

  const std::size_t n_atom =
    static_cast<std::size_t>(atom_costs.size());

  const std::size_t n_feat =
    static_cast<std::size_t>(dims[1]);

  if (static_cast<std::size_t>(dims[0]) != n_atom) {
    Rcpp::stop("feature_atom_matrix must have one row per atom.");
  }

  if (static_cast<std::size_t>(cover_by_atom.size()) != n_atom) {
    Rcpp::stop("cover_by_atom and atom_costs must have the same length.");
  }

  if (p_slot.size() != static_cast<int>(n_feat + 1)) {
    Rcpp::stop("Invalid sparse matrix column pointer slot.");
  }

  Rcpp::NumericVector targets_value = targets_list["value"];
  Rcpp::CharacterVector targets_sense = targets_list["sense"];

  if (static_cast<std::size_t>(targets_value.size()) != n_feat) {
    Rcpp::stop("Number of targets does not match number of atom feature columns.");
  }

  if (static_cast<std::size_t>(targets_sense.size()) != n_feat) {
    Rcpp::stop("Number of target senses does not match number of atom feature columns.");
  }

  // Existing rij rows are created before the objective is applied.
  // We neutralize them because this objective uses atom-based target constraints.
  std::size_t existing_rows =
    std::max(max_row_index_plus_one(ptr->_A_i), ptr->_rhs.size());

  add_neutral_rows_until(ptr, existing_rows);

  std::size_t current_ncol = n_pu;
  current_ncol = std::max(current_ncol, ptr->_obj.size());
  current_ncol = std::max(current_ncol, ptr->_lb.size());
  current_ncol = std::max(current_ncol, ptr->_ub.size());
  current_ncol = std::max(current_ncol, ptr->_vtype.size());
  current_ncol = std::max(current_ncol, ptr->_col_ids.size());

  ensure_col_vectors(ptr, current_ncol);

  // Existing PU decision variables have zero direct objective cost.
  // Cost is paid only through atom variables.
  for (std::size_t k = 0; k < current_ncol; ++k) {
    ptr->_obj[k] = 0.0;
  }

  // Add atom activation variables z_a.
  const std::size_t first_atom_col = current_ncol;

  for (std::size_t a = 0; a < n_atom; ++a) {
    ptr->_obj.push_back(atom_costs[a]);
    ptr->_lb.push_back(0.0);
    ptr->_ub.push_back(1.0);

    // Continuous is sufficient because:
    // if any covering PU is selected, constraints force z_a >= 1;
    // with ub = 1, z_a = 1.
    // if no covering PU is selected, the upper link forces z_a = 0.
    ptr->_vtype.push_back("C");
    ptr->_col_ids.push_back("ms_atom");
  }

  std::size_t row = ptr->_rhs.size();

  // Atom activation constraints:
  //
  // z_a <= sum_i x_i for PUs that cover atom a
  // x_i <= z_a for each PU i that covers atom a
  //
  // Together these make z_a active if any covering PU is selected.
  for (std::size_t a = 0; a < n_atom; ++a) {
    Rcpp::IntegerVector cover = cover_by_atom[a];
    const std::size_t z_col = first_atom_col + a;

    // z_a - sum_i x_i <= 0
    ptr->_A_i.push_back(row);
    ptr->_A_j.push_back(z_col);
    ptr->_A_x.push_back(1.0);

    for (int k = 0; k < cover.size(); ++k) {
      int i1 = cover[k];

      if (Rcpp::IntegerVector::is_na(i1) || i1 < 1) {
        Rcpp::stop("Invalid planning-unit index in cover_by_atom.");
      }

      std::size_t i = static_cast<std::size_t>(i1 - 1);

      if (i >= n_pu) {
        Rcpp::stop("cover_by_atom contains out-of-range planning-unit index.");
      }

      ptr->_A_i.push_back(row);
      ptr->_A_j.push_back(i);
      ptr->_A_x.push_back(-1.0);
    }

    ptr->_sense.push_back("<=");
    ptr->_rhs.push_back(0.0);
    ptr->_row_ids.push_back("ms_atom_upper");
    ++row;

    // x_i - z_a <= 0
    for (int k = 0; k < cover.size(); ++k) {
      int i1 = cover[k];

      if (Rcpp::IntegerVector::is_na(i1) || i1 < 1) {
        Rcpp::stop("Invalid planning-unit index in cover_by_atom.");
      }

      std::size_t i = static_cast<std::size_t>(i1 - 1);

      if (i >= n_pu) {
        Rcpp::stop("cover_by_atom contains out-of-range planning-unit index.");
      }

      ptr->_A_i.push_back(row);
      ptr->_A_j.push_back(i);
      ptr->_A_x.push_back(1.0);

      ptr->_A_i.push_back(row);
      ptr->_A_j.push_back(z_col);
      ptr->_A_x.push_back(-1.0);

      ptr->_sense.push_back("<=");
      ptr->_rhs.push_back(0.0);
      ptr->_row_ids.push_back("ms_atom_lower");
      ++row;
    }
  }

  // Atom-based feature target constraints:
  //
  // sum_a q_aj z_a >= T_j
  //
  // The dgCMatrix is column-compressed:
  // column j entries are k = p[j], ..., p[j + 1] - 1
  for (std::size_t j = 0; j < n_feat; ++j) {
    const int start = p_slot[j];
    const int end = p_slot[j + 1];

    for (int k = start; k < end; ++k) {
      const int a_int = i_slot[k];

      if (a_int < 0) {
        Rcpp::stop("Invalid negative row index in sparse matrix.");
      }

      const std::size_t a = static_cast<std::size_t>(a_int);

      if (a >= n_atom) {
        Rcpp::stop("Sparse matrix row index exceeds number of atoms.");
      }

      const double val = x_slot[k];

      if (val == 0.0) {
        continue;
      }

      ptr->_A_i.push_back(row);
      ptr->_A_j.push_back(first_atom_col + a);
      ptr->_A_x.push_back(val);
    }

    std::string s = Rcpp::as<std::string>(targets_sense[j]);

    ptr->_sense.push_back(s);
    ptr->_rhs.push_back(targets_value[j]);
    ptr->_row_ids.push_back("ms_feature_target");
    ++row;
  }

  ptr->_modelsense = "min";
  // Preserve the first and last constraint rows when reconstructing the sparse
  // optimization matrix, including rows with no nonzero coefficients.
  if (row > 0 && !ptr->_obj.empty()) {
    ptr->_A_i.push_back(0); ptr->_A_j.push_back(0); ptr->_A_x.push_back(0.0);
    ptr->_A_i.push_back(row - 1); ptr->_A_j.push_back(0); ptr->_A_x.push_back(0.0);
  }
  Rcpp::List out = Rcpp::clone(model);
  out["modelsense"] = ptr->_modelsense;
  out["number_of_features"] = ptr->_number_of_features;
  out["number_of_planning_units"] = ptr->_number_of_planning_units;
  out["number_of_zones"] = ptr->_number_of_zones;
  out["A_i"] = Rcpp::IntegerVector(ptr->_A_i.begin(), ptr->_A_i.end());
  out["A_j"] = Rcpp::IntegerVector(ptr->_A_j.begin(), ptr->_A_j.end());
  out["A_x"] = ptr->_A_x;
  out["obj"] = ptr->_obj;
  out["lb"] = ptr->_lb;
  out["ub"] = ptr->_ub;
  out["rhs"] = ptr->_rhs;
  out["sense"] = ptr->_sense;
  out["vtype"] = ptr->_vtype;
  out["row_ids"] = ptr->_row_ids;
  out["col_ids"] = ptr->_col_ids;
  out["compressed_formulation"] = ptr->_compressed_formulation;
  return out;
}
