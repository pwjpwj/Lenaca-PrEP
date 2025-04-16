#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List HIVSIM_cpp(int Age, int NoSem, int NumHSH, int NumHSHVIH,
                double ACD4, double I8CD4, double I7CD4, double Rest,
                double G, double Woff, double W, NumericVector m,
                double mVIH235, double mVIH35, double mAIDS,
                double mTr, double Cvh, double Ch, double Cl, double Cvl,
                double nvh, double nh, double nl, double nvl,
                double B, double I12, double I23, double I34,
                double I45, double I55, double I56, double I67,
                double I78, double I89, NumericVector PTr, int years,
                NumericVector PropTR, NumericVector DistEdad, double RiskRed) {
  
  List Q(years);
  
  NumericVector Compartments(NoSem * Age * 10 * 4 * 2, 0.0);
  IntegerVector dim_compartments = {NoSem, Age, 10, 4, 2};
  Compartments.attr("dim") = dim_compartments;
  
  NumericMatrix Tr(NoSem, Age, 0.0);
  NumericMatrix Dead(NoSem, Age, 0.0);
  NumericMatrix Lvh_mat(NoSem, Age, 0.0);
  NumericMatrix Lh_mat(NoSem, Age, 0.0);
  NumericMatrix Ll_mat(NoSem, Age, 0.0);
  NumericMatrix Lvl_mat(NoSem, Age, 0.0);
  NumericMatrix CalAIDS(NoSem, Age, 0.0);
  NumericMatrix NewHIV(NoSem, Age, 0.0);
  NumericMatrix NewTreatments(NoSem, Age, 0.0);
  
  for (int j = 0; j < Age; ++j) {
    for (int k = 0; k < 4; ++k) {
      std::vector<double> RG = {static_cast<double>(nvh), static_cast<double>(nh), 
                                static_cast<double>(nl), static_cast<double>(nvl)};
      
      if (j >= DistEdad.size()) continue;
      
      Compartments.at(j + k * Age * 10 * 2) = RG.at(k) * NumHSH * DistEdad.at(j);
      
      for (int stage = 1; stage <= 6; ++stage) {
        Compartments.at(j + k * Age * 10 * 2 + stage * Age * 10 * 4 * 2) =
          NumHSHVIH * RG.at(k) * Rest / 6.0 * DistEdad.at(j);
      }
      
      for (int stage = 7; stage <= 8; ++stage) {
        Compartments.at(j + k * Age * 10 * 2 + stage * Age * 10 * 4 * 2) =
          NumHSHVIH * I7CD4 * RG.at(k) * DistEdad.at(j);
      }
      
      Compartments.at(j + k * Age * 10 * 2 + 9 * Age * 10 * 4 * 2) =
        NumHSHVIH * RG.at(k) * ACD4 * DistEdad.at(j);
    }
    
    Tr(0, j) = NumHSHVIH * PropTR.at(j);
  }
  
  for (int t = 0; t < years; ++t) {
    for (int j = 0; j < Age; ++j) {
      if (t >= 1 && j < 84) {
        for (int k = 0; k < 4; ++k) {
          for (int prep = 0; prep < 2; ++prep) {
            for (int stage = 0; stage < 10; ++stage) {
              int index_t1 = j + k * Age * 10 * 2 + stage * Age * 10 * 4 * 2 + prep * Age * 10 * 4 * 2 * NoSem;
              int index_t0 = (j + 1) + k * Age * 10 * 2 + stage * Age * 10 * 4 * 2 + prep * Age * 10 * 4 * 2 * NoSem + 51 * Age * 10 * 4 * 2;
              Compartments.at(index_t1) = Compartments.at(index_t0);
            }
          }
        }
        Tr(0, j + 1) = Tr(51, j);
      }
    
    
    IntegerVector week_col(NoSem);
    for (int week = 0; week < NoSem; ++week) {
      week_col[week] = week + 1;
    }
    
    NumericMatrix Compartments_matrix(NoSem, 10 * 4 * 2);
    for (int week = 0; week < NoSem; ++week) {
      for (int stage = 0; stage < 10; ++stage) {
        for (int risk = 0; risk < 4; ++risk) {
          for (int prep = 0; prep < 2; ++prep) {
            int index = week * Age * 10 * 4 * 2 + j + stage * Age * 10 * 4 * 2 + risk * Age * 10 * 2 + prep * Age * 10;
            Compartments_matrix(week, stage * 4 * 2 + risk * 2 + prep) = Compartments.at(index);
          }
        }
      }
    }
    
    DataFrame P = DataFrame::create(
      Named("Week") = week,
      Named("Compartments") = Compartments,
      Named("Tr") = Tr,
      Named("Dead") = Dead,
      Named("CalAIDS") = CalAIDS,
      Named("Lvh_mat") = Lvh_mat,
      Named("Lh_mat") = Lh_mat,
      Named("Ll_mat") = Ll_mat,
      Named("Lvl_mat") = Lvl_mat,
      Named("NewHIV") = NewHIV,
      Named("NewTreatments") = NewTreatments
    );
    
    Q[t] = P;
  }
  }
  return Q;
}
