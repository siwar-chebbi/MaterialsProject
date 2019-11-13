#!/usr/bin/env python

# Test script for gbml elasticity (bulk and shear moduli) predictions

from gbml import elasticity

API_KEY = "eDCEK5m9WVjmajp7e8af"

mpID = "mp-10003"
(k_value, g_value, caveat_str) = elasticity.predict_k_g(mpID, api_key=API_KEY)
print (k_value, g_value, caveat_str)

mpID_list = ["mp-10003","mp-10010","mp-10015","mp-10021","mp-26","mp-10018","mp-19306"]
(matid_list, k_list, g_list, caveat_list) = elasticity.predict_k_g_list(mpID_list, api_key=API_KEY)
print (matid_list, k_list, g_list, caveat_list)
