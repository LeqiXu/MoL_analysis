for file in {1..200}
do
if [[ ! -e "${file}.txt" ]]; then
echo ${file}
fi
done

# only consider n = 2000, 5000, 10000, 20000, 50000, 100000 p = 1000 w_0.15_h_0.4 in norm

# dsq-MoL_mix_n_2e+05_p_20000_w_0.375_h_0.8_N_10_mu_0.2
# dsq-MoL_mix_n_5e+05_p_10000_w_0.375_h_0.8_N_10_mu_0.2
# dsq-MoL_mix_n_5e+05_p_10000_w_0.0625_h_0.8_N_10_mu_0.2
# dsq-MoL_mix_n_5e+05_p_10000_w_0.5_h_0.4_N_10_mu_0.2
# dsq-MoL_mix_n_3e+05_p_5000_w_0.15_h_0.4_N_10_mu_0.2
# dsq-MoL_mix_n_3e+05_p_5000_w_0.5_h_0.4_N_10_mu_0.2
# dsq-MoL_mix_n_3e+05_p_5000_w_0.0625_h_0.8_N_10_mu_0.2
# dsq-MoL_mix_n_5e+05_p_10000_w_0.15_h_0.4_N_10_mu_0.2
# dsq-MoL_mix_n_5e+05_p_5000_w_0.375_h_0.8_N_10_mu_0.2
# dsq-MoL_mix_n_5e+05_p_5000_w_0.0625_h_0.8_N_10_mu_0.2
# dsq-MoL_mix_n_5e+05_p_5000_w_0.5_h_0.4_N_10_mu_0.2
# dsq-MoL_mix_n_5e+05_p_5000_w_0.15_h_0.4_N_10_mu_0.2
# dsq-MoL_mix_n_3e+05_p_20000_w_0.375_h_0.8_N_10_mu_0.2
# dsq-MoL_mix_n_3e+05_p_20000_w_0.0625_h_0.8_N_10_mu_0.2
# dsq-MoL_mix_n_3e+05_p_20000_w_0.5_h_0.4_N_10_mu_0.2
# dsq-MoL_mix_n_3e+05_p_20000_w_0.15_h_0.4_N_10_mu_0.2
# dsq-MoL_mix_n_3e+05_p_10000_w_0.375_h_0.8_N_10_mu_0.2
# dsq-MoL_mix_n_3e+05_p_10000_w_0.0625_h_0.8_N_10_mu_0.2
# dsq-MoL_mix_n_3e+05_p_10000_w_0.5_h_0.4_N_10_mu_0.2
# dsq-MoL_mix_n_3e+05_p_10000_w_0.15_h_0.4_N_10_mu_0.2
# dsq-MoL_mix_n_3e+05_p_5000_w_0.375_h_0.8_N_10_mu_0.2


vim MoL_rest.sh
Rscript --vanilla General_MoL.R mix 50000 20000 0.5 0.4 1
Rscript --vanilla General_MoL.R mix 50000 20000 0.5 0.4 63
Rscript --vanilla General_MoL.R mix 50000 20000 0.5 0.4 64
Rscript --vanilla General_MoL.R mix 50000 20000 0.5 0.4 74
Rscript --vanilla General_MoL.R mix 50000 20000 0.5 0.4 75
Rscript --vanilla General_MoL.R mix 50000 20000 0.5 0.4 79
Rscript --vanilla General_MoL.R mix 50000 20000 0.5 0.4 162
Rscript --vanilla General_MoL.R mix 50000 20000 0.5 0.4 163
Rscript --vanilla General_MoL.R mix 50000 20000 0.5 0.4 173
Rscript --vanilla General_MoL.R mix 50000 20000 0.5 0.4 174
Rscript --vanilla General_MoL.R mix 50000 20000 0.5 0.4 188
Rscript --vanilla General_MoL.R mix 50000 20000 0.5 0.4 198
Rscript --vanilla General_MoL.R mix 50000 20000 0.5 0.4 199
Rscript --vanilla General_MoL.R mix 50000 10000 0.375 0.8 156
Rscript --vanilla General_MoL.R mix 50000 10000 0.375 0.8 171
Rscript --vanilla General_MoL.R mix 50000 20000 0.0625 0.8 39
Rscript --vanilla General_MoL.R mix 50000 20000 0.0625 0.8 40
Rscript --vanilla General_MoL.R mix 50000 20000 0.0625 0.8 41
Rscript --vanilla General_MoL.R mix 50000 20000 0.0625 0.8 42
Rscript --vanilla General_MoL.R mix 50000 20000 0.0625 0.8 74
Rscript --vanilla General_MoL.R mix 50000 20000 0.0625 0.8 75
Rscript --vanilla General_MoL.R mix 50000 20000 0.0625 0.8 145
Rscript --vanilla General_MoL.R mix 50000 20000 0.0625 0.8 148
Rscript --vanilla General_MoL.R mix 50000 20000 0.0625 0.8 150
Rscript --vanilla General_MoL.R mix 50000 20000 0.0625 0.8 151
Rscript --vanilla General_MoL.R mix 50000 20000 0.0625 0.8 152
Rscript --vanilla General_MoL.R mix 50000 20000 0.0625 0.8 161
Rscript --vanilla General_MoL.R mix 50000 20000 0.0625 0.8 162
Rscript --vanilla General_MoL.R mix 50000 20000 0.0625 0.8 163
Rscript --vanilla General_MoL.R mix 50000 20000 0.0625 0.8 164


Rscript --vanilla General_MoL.R mix 1e+05 10000 0.375 0.8 81
Rscript --vanilla General_MoL.R mix 1e+05 10000 0.375 0.8 83
Rscript --vanilla General_MoL.R mix 1e+05 10000 0.375 0.8 110
Rscript --vanilla General_MoL.R mix 1e+05 10000 0.375 0.8 112
Rscript --vanilla General_MoL.R mix 1e+05 10000 0.375 0.8 113
Rscript --vanilla General_MoL.R mix 1e+05 10000 0.375 0.8 131
Rscript --vanilla General_MoL.R mix 1e+05 10000 0.375 0.8 137
Rscript --vanilla General_MoL.R mix 1e+05 10000 0.5 0.4 40
Rscript --vanilla General_MoL.R mix 1e+05 10000 0.5 0.4 41
Rscript --vanilla General_MoL.R mix 1e+05 10000 0.5 0.4 44
Rscript --vanilla General_MoL.R mix 1e+05 10000 0.5 0.4 105
Rscript --vanilla General_MoL.R mix 1e+05 10000 0.5 0.4 110
Rscript --vanilla General_MoL.R mix 1e+05 10000 0.5 0.4 113
Rscript --vanilla General_MoL.R mix 1e+05 10000 0.5 0.4 116
Rscript --vanilla General_MoL.R mix 2e+05 5000 0.0625 0.8 49
Rscript --vanilla General_MoL.R mix 2e+05 5000 0.0625 0.8 52
Rscript --vanilla General_MoL.R mix 2e+05 5000 0.0625 0.8 155
Rscript --vanilla General_MoL.R mix 2e+05 5000 0.5 0.4 101
Rscript --vanilla General_MoL.R mix 2e+05 5000 0.5 0.4 105
Rscript --vanilla General_MoL.R mix 2e+05 5000 0.5 0.4 114
Rscript --vanilla General_MoL.R mix 2e+05 5000 0.5 0.4 140
Rscript --vanilla General_MoL.R mix 2e+05 5000 0.15 0.4 4
Rscript --vanilla General_MoL.R mix 3e+05 2000 0.375 0.8 1
Rscript --vanilla General_MoL.R mix 3e+05 2000 0.375 0.8 95
Rscript --vanilla General_MoL.R mix 3e+05 2000 0.375 0.8 179
Rscript --vanilla General_MoL.R mix 3e+05 2000 0.375 0.8 180
Rscript --vanilla General_MoL.R mix 3e+05 2000 0.375 0.8 181
Rscript --vanilla General_MoL.R mix 3e+05 2000 0.375 0.8 198
Rscript --vanilla General_MoL.R mix 3e+05 2000 0.0625 0.8 123
Rscript --vanilla General_MoL.R mix 3e+05 2000 0.5 0.4 22
Rscript --vanilla General_MoL.R mix 3e+05 2000 0.5 0.4 23
Rscript --vanilla General_MoL.R mix 3e+05 2000 0.15 0.4 1
Rscript --vanilla General_MoL.R mix 3e+05 2000 0.15 0.4 2
Rscript --vanilla General_MoL.R mix 3e+05 2000 0.15 0.4 3
Rscript --vanilla General_MoL.R mix 3e+05 2000 0.15 0.4 4
Rscript --vanilla General_MoL.R mix 3e+05 2000 0.15 0.4 5
Rscript --vanilla General_MoL.R mix 3e+05 1000 0.375 0.8 64
Rscript --vanilla General_MoL.R mix 3e+05 1000 0.375 0.8 65
Rscript --vanilla General_MoL.R mix 3e+05 1000 0.375 0.8 82
Rscript --vanilla General_MoL.R mix 3e+05 1000 0.375 0.8 83
Rscript --vanilla General_MoL.R mix 3e+05 1000 0.375 0.8 84
Rscript --vanilla General_MoL.R mix 3e+05 1000 0.375 0.8 85
Rscript --vanilla General_MoL.R mix 3e+05 1000 0.375 0.8 86
Rscript --vanilla General_MoL.R mix 3e+05 1000 0.375 0.8 87
Rscript --vanilla General_MoL.R mix 3e+05 1000 0.375 0.8 100
Rscript --vanilla General_MoL.R mix 3e+05 1000 0.375 0.8 101
Rscript --vanilla General_MoL.R mix 3e+05 1000 0.375 0.8 102
Rscript --vanilla General_MoL.R mix 3e+05 1000 0.375 0.8 165
Rscript --vanilla General_MoL.R mix 3e+05 1000 0.375 0.8 166
Rscript --vanilla General_MoL.R mix 3e+05 1000 0.375 0.8 167
Rscript --vanilla General_MoL.R mix 3e+05 1000 0.375 0.8 168
Rscript --vanilla General_MoL.R mix 3e+05 1000 0.375 0.8 169
Rscript --vanilla General_MoL.R mix 3e+05 1000 0.375 0.8 170
Rscript --vanilla General_MoL.R mix 3e+05 1000 0.375 0.8 171
Rscript --vanilla General_MoL.R mix 3e+05 1000 0.375 0.8 172
Rscript --vanilla General_MoL.R mix 3e+05 1000 0.375 0.8 182
Rscript --vanilla General_MoL.R mix 3e+05 1000 0.375 0.8 183
Rscript --vanilla General_MoL.R mix 3e+05 1000 0.375 0.8 184
Rscript --vanilla General_MoL.R mix 3e+05 1000 0.375 0.8 185
Rscript --vanilla General_MoL.R mix 3e+05 10000 0.15 0.4 57
Rscript --vanilla General_MoL.R mix 5e+05 2000 0.375 0.8 24
Rscript --vanilla General_MoL.R mix 5e+05 2000 0.375 0.8 52
Rscript --vanilla General_MoL.R mix 5e+05 2000 0.0625 0.8 1
Rscript --vanilla General_MoL.R mix 5e+05 2000 0.0625 0.8 14
Rscript --vanilla General_MoL.R mix 5e+05 2000 0.5 0.4 7
Rscript --vanilla General_MoL.R mix 5e+05 2000 0.5 0.4 187
Rscript --vanilla General_MoL.R norm 3e+05 10000 0.0625 0.8 55
Rscript --vanilla General_MoL.R norm 3e+05 10000 0.0625 0.8 125
Rscript --vanilla General_MoL.R norm 3e+05 10000 0.0625 0.8 161
Rscript --vanilla General_MoL.R norm 3e+05 10000 0.0625 0.8 189
Rscript --vanilla General_MoL.R norm 3e+05 10000 0.0625 0.8 190
Rscript --vanilla General_MoL.R norm 3e+05 10000 0.0625 0.8 191
Rscript --vanilla General_MoL.R norm 3e+05 10000 0.0625 0.8 192
Rscript --vanilla General_MoL.R norm 3e+05 10000 0.0625 0.8 193
Rscript --vanilla General_MoL.R norm 3e+05 10000 0.0625 0.8 194
Rscript --vanilla General_MoL.R norm 3e+05 10000 0.0625 0.8 195
Rscript --vanilla General_MoL.R norm 3e+05 10000 0.0625 0.8 196
Rscript --vanilla General_MoL.R norm 3e+05 10000 0.5 0.4 17
Rscript --vanilla General_MoL.R norm 3e+05 10000 0.5 0.4 18
Rscript --vanilla General_MoL.R norm 3e+05 10000 0.5 0.4 20
Rscript --vanilla General_MoL.R norm 3e+05 10000 0.5 0.4 23
Rscript --vanilla General_MoL.R norm 3e+05 10000 0.5 0.4 27
Rscript --vanilla General_MoL.R norm 3e+05 10000 0.5 0.4 28
Rscript --vanilla General_MoL.R norm 3e+05 10000 0.5 0.4 29
Rscript --vanilla General_MoL.R norm 3e+05 10000 0.5 0.4 34
Rscript --vanilla General_MoL.R norm 3e+05 10000 0.5 0.4 42
Rscript --vanilla General_MoL.R norm 3e+05 10000 0.5 0.4 43
Rscript --vanilla General_MoL.R norm 3e+05 10000 0.5 0.4 45
Rscript --vanilla General_MoL.R norm 3e+05 10000 0.5 0.4 46
Rscript --vanilla General_MoL.R norm 3e+05 10000 0.5 0.4 47
Rscript --vanilla General_MoL.R norm 3e+05 10000 0.5 0.4 48
Rscript --vanilla General_MoL.R norm 3e+05 10000 0.375 0.8 6
Rscript --vanilla General_MoL.R norm 3e+05 10000 0.375 0.8 35
Rscript --vanilla General_MoL.R norm 3e+05 10000 0.375 0.8 36
Rscript --vanilla General_MoL.R norm 3e+05 10000 0.375 0.8 37
Rscript --vanilla General_MoL.R norm 3e+05 10000 0.375 0.8 48
Rscript --vanilla General_MoL.R norm 3e+05 10000 0.375 0.8 56
Rscript --vanilla General_MoL.R norm 3e+05 10000 0.375 0.8 59
Rscript --vanilla General_MoL.R norm 3e+05 10000 0.375 0.8 93
Rscript --vanilla General_MoL.R norm 3e+05 10000 0.375 0.8 100
Rscript --vanilla General_MoL.R norm 3e+05 10000 0.375 0.8 101
Rscript --vanilla General_MoL.R norm 3e+05 10000 0.375 0.8 102
Rscript --vanilla General_MoL.R norm 3e+05 10000 0.375 0.8 188
Rscript --vanilla General_MoL.R norm 3e+05 10000 0.375 0.8 189

Rscript --vanilla General_MoL.R mix 1e+05 20000 0.15 0.4 97
Rscript --vanilla General_MoL.R mix 1e+05 20000 0.15 0.4 98
Rscript --vanilla General_MoL.R mix 1e+05 20000 0.15 0.4 99
Rscript --vanilla General_MoL.R mix 1e+05 20000 0.15 0.4 100
Rscript --vanilla General_MoL.R mix 1e+05 20000 0.15 0.4 103
Rscript --vanilla General_MoL.R mix 1e+05 20000 0.15 0.4 104
Rscript --vanilla General_MoL.R mix 1e+05 20000 0.15 0.4 105
Rscript --vanilla General_MoL.R mix 1e+05 20000 0.15 0.4 106
Rscript --vanilla General_MoL.R mix 1e+05 20000 0.15 0.4 107
Rscript --vanilla General_MoL.R mix 1e+05 20000 0.15 0.4 147
Rscript --vanilla General_MoL.R mix 1e+05 20000 0.15 0.4 149
Rscript --vanilla General_MoL.R mix 1e+05 20000 0.15 0.4 151
Rscript --vanilla General_MoL.R mix 1e+05 20000 0.15 0.4 153
Rscript --vanilla General_MoL.R mix 1e+05 20000 0.15 0.4 158
Rscript --vanilla General_MoL.R mix 1e+05 20000 0.15 0.4 159
Rscript --vanilla General_MoL.R mix 1e+05 20000 0.15 0.4 162
Rscript --vanilla General_MoL.R mix 1e+05 20000 0.15 0.4 193
Rscript --vanilla General_MoL.R mix 1e+05 20000 0.15 0.4 195
Rscript --vanilla General_MoL.R mix 1e+05 20000 0.15 0.4 196
Rscript --vanilla General_MoL.R mix 1e+05 20000 0.15 0.4 197
Rscript --vanilla General_MoL.R mix 1e+05 20000 0.15 0.4 200

Rscript --vanilla General_MoL.R norm 5e+05 10000 0.5 0.4 155
Rscript --vanilla General_MoL.R norm 5e+05 10000 0.5 0.4 156
Rscript --vanilla General_MoL.R norm 5e+05 10000 0.5 0.4 157
Rscript --vanilla General_MoL.R norm 5e+05 10000 0.5 0.4 158
Rscript --vanilla General_MoL.R norm 5e+05 10000 0.15 0.4 132
Rscript --vanilla General_MoL.R norm 5e+05 10000 0.15 0.4 138
Rscript --vanilla General_MoL.R norm 5e+05 10000 0.15 0.4 143
Rscript --vanilla General_MoL.R norm 5e+05 10000 0.15 0.4 147
Rscript --vanilla General_MoL.R norm 5e+05 10000 0.15 0.4 156
Rscript --vanilla General_MoL.R norm 5e+05 10000 0.15 0.4 168
Rscript --vanilla General_MoL.R norm 3e+05 20000 0.0625 0.8 26
Rscript --vanilla General_MoL.R norm 3e+05 20000 0.0625 0.8 99
Rscript --vanilla General_MoL.R norm 3e+05 20000 0.0625 0.8 107
Rscript --vanilla General_MoL.R norm 3e+05 20000 0.5 0.4 158
Rscript --vanilla General_MoL.R norm 3e+05 20000 0.5 0.4 180
Rscript --vanilla General_MoL.R norm 2e+05 20000 0.15 0.4 81
Rscript --vanilla General_MoL.R norm 2e+05 20000 0.15 0.4 100
Rscript --vanilla General_MoL.R mix 2e+05 20000 0.15 0.4 82
Rscript --vanilla General_MoL.R mix 2e+05 20000 0.15 0.4 88
Rscript --vanilla General_MoL.R mix 2e+05 20000 0.15 0.4 93
Rscript --vanilla General_MoL.R mix 2e+05 20000 0.15 0.4 95
Rscript --vanilla General_MoL.R mix 2e+05 20000 0.15 0.4 96
Rscript --vanilla General_MoL.R mix 2e+05 20000 0.15 0.4 98
Rscript --vanilla General_MoL.R mix 2e+05 20000 0.15 0.4 100
Rscript --vanilla General_MoL.R mix 2e+05 20000 0.15 0.4 114
Rscript --vanilla General_MoL.R mix 2e+05 20000 0.15 0.4 130
Rscript --vanilla General_MoL.R mix 2e+05 20000 0.15 0.4 132
Rscript --vanilla General_MoL.R mix 2e+05 20000 0.15 0.4 178
Rscript --vanilla General_MoL.R mix 2e+05 20000 0.15 0.4 179
Rscript --vanilla General_MoL.R mix 2e+05 20000 0.15 0.4 185
Rscript --vanilla General_MoL.R mix 2e+05 20000 0.15 0.4 186
Rscript --vanilla General_MoL.R mix 2e+05 20000 0.15 0.4 199
Rscript --vanilla General_MoL.R mix 2e+05 20000 0.5 0.4 4
Rscript --vanilla General_MoL.R mix 2e+05 20000 0.5 0.4 8
Rscript --vanilla General_MoL.R mix 2e+05 20000 0.5 0.4 9
Rscript --vanilla General_MoL.R mix 2e+05 20000 0.5 0.4 10
Rscript --vanilla General_MoL.R mix 2e+05 20000 0.5 0.4 11
Rscript --vanilla General_MoL.R mix 2e+05 20000 0.5 0.4 12
Rscript --vanilla General_MoL.R mix 2e+05 20000 0.5 0.4 19
Rscript --vanilla General_MoL.R mix 2e+05 20000 0.5 0.4 27
Rscript --vanilla General_MoL.R mix 2e+05 20000 0.5 0.4 43
Rscript --vanilla General_MoL.R mix 2e+05 20000 0.5 0.4 44
Rscript --vanilla General_MoL.R mix 2e+05 20000 0.5 0.4 67
Rscript --vanilla General_MoL.R mix 2e+05 20000 0.5 0.4 73
Rscript --vanilla General_MoL.R mix 2e+05 20000 0.5 0.4 74
Rscript --vanilla General_MoL.R mix 2e+05 20000 0.5 0.4 75
Rscript --vanilla General_MoL.R mix 2e+05 20000 0.5 0.4 77
Rscript --vanilla General_MoL.R mix 2e+05 20000 0.5 0.4 78
Rscript --vanilla General_MoL.R mix 2e+05 20000 0.5 0.4 90
Rscript --vanilla General_MoL.R mix 2e+05 20000 0.5 0.4 100
Rscript --vanilla General_MoL.R mix 2e+05 20000 0.5 0.4 105
Rscript --vanilla General_MoL.R mix 2e+05 20000 0.5 0.4 107
Rscript --vanilla General_MoL.R mix 2e+05 20000 0.5 0.4 113
Rscript --vanilla General_MoL.R mix 2e+05 20000 0.5 0.4 167
Rscript --vanilla General_MoL.R mix 2e+05 20000 0.5 0.4 168
Rscript --vanilla General_MoL.R mix 2e+05 20000 0.5 0.4 172
Rscript --vanilla General_MoL.R mix 2e+05 20000 0.5 0.4 174
Rscript --vanilla General_MoL.R mix 2e+05 20000 0.5 0.4 176
Rscript --vanilla General_MoL.R mix 2e+05 20000 0.5 0.4 179
Rscript --vanilla General_MoL.R mix 2e+05 20000 0.5 0.4 180
Rscript --vanilla General_MoL.R mix 2e+05 20000 0.5 0.4 185
Rscript --vanilla General_MoL.R mix 2e+05 20000 0.5 0.4 191
Rscript --vanilla General_MoL.R mix 2e+05 20000 0.5 0.4 197
Rscript --vanilla General_MoL.R mix 2e+05 20000 0.5 0.4 200

vim PQL_norm_rest.sh
Rscript --vanilla General_PQLseq.R norm 20000 1000 0.15 0.4 109
Rscript --vanilla General_PQLseq.R norm 20000 1000 0.15 0.4 174

vim PQL_norm_rest.submit
module load GCC;
partition=scavenge;
module load dSQ
dSQ --jobfile PQL_norm_rest.sh -p ${partition} -n 1 --constraint avx2 --mem-per-cpu=32g -t 24:00:00 --mail-type=ALL --batch-file PQL_mix_rest.pbs
sbatch PQL_norm_rest.pbs

vim PQL_mix_rest.sh
Rscript --vanilla General_PQLseq.R mix 1000 2000 0.15 0.4 148
Rscript --vanilla General_PQLseq.R mix 5000 500 0.375 0.8 146
Rscript --vanilla General_PQLseq.R mix 5000 500 0.375 0.8 151
Rscript --vanilla General_PQLseq.R mix 5000 500 0.375 0.8 159
Rscript --vanilla General_PQLseq.R mix 5000 500 0.375 0.8 162
Rscript --vanilla General_PQLseq.R mix 5000 500 0.375 0.8 163
Rscript --vanilla General_PQLseq.R mix 5000 500 0.375 0.8 166
Rscript --vanilla General_PQLseq.R mix 5000 500 0.375 0.8 169
Rscript --vanilla General_PQLseq.R mix 5000 500 0.375 0.8 170
Rscript --vanilla General_PQLseq.R mix 5000 500 0.375 0.8 178
Rscript --vanilla General_PQLseq.R mix 5000 500 0.375 0.8 194
Rscript --vanilla General_PQLseq.R mix 5000 500 0.15 0.4 107
Rscript --vanilla General_PQLseq.R mix 5000 500 0.15 0.4 111
Rscript --vanilla General_PQLseq.R mix 5000 500 0.15 0.4 115
Rscript --vanilla General_PQLseq.R mix 5000 500 0.15 0.4 117
Rscript --vanilla General_PQLseq.R mix 5000 500 0.15 0.4 124
Rscript --vanilla General_PQLseq.R mix 5000 500 0.15 0.4 127
Rscript --vanilla General_PQLseq.R mix 5000 500 0.15 0.4 134
Rscript --vanilla General_PQLseq.R mix 5000 500 0.15 0.4 144
Rscript --vanilla General_PQLseq.R mix 5000 500 0.15 0.4 148
Rscript --vanilla General_PQLseq.R mix 5000 500 0.15 0.4 149
Rscript --vanilla General_PQLseq.R mix 5000 500 0.15 0.4 160
Rscript --vanilla General_PQLseq.R mix 5000 500 0.15 0.4 165
Rscript --vanilla General_PQLseq.R mix 5000 500 0.15 0.4 183
Rscript --vanilla General_PQLseq.R mix 2000 500 0.15 0.4 118
Rscript --vanilla General_PQLseq.R mix 2000 500 0.15 0.4 119
Rscript --vanilla General_PQLseq.R mix 2000 500 0.15 0.4 120
Rscript --vanilla General_PQLseq.R mix 2000 500 0.15 0.4 129
Rscript --vanilla General_PQLseq.R mix 2000 500 0.15 0.4 131
Rscript --vanilla General_PQLseq.R mix 2000 500 0.15 0.4 135
Rscript --vanilla General_PQLseq.R mix 2000 500 0.15 0.4 136
Rscript --vanilla General_PQLseq.R mix 2000 500 0.15 0.4 143
Rscript --vanilla General_PQLseq.R mix 2000 500 0.15 0.4 144
Rscript --vanilla General_PQLseq.R mix 2000 500 0.15 0.4 146

vim PQL_mix_rest.submit
module load GCC;
partition=scavenge;
module load dSQ
dSQ --jobfile PQL_mix_rest.sh -p ${partition} -n 1 --constraint avx2 --mem-per-cpu=32g -t 24:00:00 --mail-type=ALL --batch-file PQL_mix_rest.pbs
sbatch PQL_mix_rest.pbs

vim submit_one_PQL.sh
#!/bin/bash
#SBATCH --partition=bigmem
#SBATCH --mem=1000G
#SBATCH --cpus-per-task=5
#SBATCH --ntasks=1 --nodes=1
#SBATCH --time=3-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name=PQL_n_20000_p_1000_w_0.15_h_0.4
#SBATCH --output=out_PQL_one.txt

module load R

## Method
Rscript --vanilla General_PQLseq.R norm 20000 1000 0.15 0.4 1