/* EventMatching
 * event_matching.cpp
 * This program contains an adaptation of the stereo_matching_bp.m algorithm 
 * (obtained from https://github.com/harryxz/EMP/releases/tag/Frontiers).
 * It is based on the paper by Xie et al. 
 * (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5633728/).
 *  Created on: 24 May 2018
 *      Authors: Nico Hertel, Alexander Matthies, Seth Siriya
 */

// set constants
#define EVENT_SIZE 6
#define TS_COL 0
#define X_COL 1
#define Y_COL 2
#define P_COL 3
#define C_COL 4
#define GT_COL 5
#define LEFT_CAM 0
#define RIGHT_CAM 1
#define OLD_MRF_THRESHOLD 10000

#include <iostream>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <armadillo>
#include <ctime>

using namespace std;
using namespace arma;

arma::Mat<double> stereo_matching_bp(arma::Mat<double> &TD, bool groundtruth_bool);

int main() {

	string csv_filename;
	string out_filename;
	cout << "Enter csv filename: ";
	cin >> csv_filename;
	cout << "Enter output filename: ";
	cin >> out_filename;

	// Load event data into matrix
	arma::Mat<double> TD;
	TD.load(csv_filename, arma::csv_ascii);

	// Apply event matching
	time_t start = time(0);
	arma::Mat<double> stereo_TD = stereo_matching_bp(TD, true);
	time_t end = time(0);

	// Output to csv file
	stereo_TD.save(out_filename, csv_ascii);

	// Show time
	double time = difftime(end, start);
	cout << "elapsed time (s): " << time << endl;

	return 0;
}

arma::Mat<double> stereo_matching_bp(arma::Mat<double> &TD, bool groundtruth_bool) {

	int disparity;

	// parameters for EMP
	int width = 240;
	int height = 180;
	double temporal_var_Thre = 3000;
	double special_var_Thre = 3;
	int wradius = 1;
	int LABEL = 50;
	double THETA = 1;
	double DISC_K = 1.7;
	double DATE_TRUNC = 5;
	double BP_ITERATIONS = 1;
	double LAMBDA = 1;
	int WRADIUS_TS = 2;

	double old_pixel_threshold = 20000;

	// last spike map store the latest ts and p
	arma::Cube<double> last_spike_time_right(height, width, 2);
	last_spike_time_right.fill(datum::nan);

	// define the struct of the MRF field and set messages to zero
	arma::Cube<double> msg_up(height, width, LABEL, fill::zeros);
	arma::Cube<double> msg_down(height, width, LABEL, fill::zeros);
	arma::Cube<double> msg_left(height, width, LABEL, fill::zeros);
	arma::Cube<double> msg_right(height, width, LABEL, fill::zeros);
	arma::Cube<double> msg_data(height, width, LABEL, fill::zeros);
	arma::Cube<double> msg_ts(height, width, LABEL, fill::zeros);
	msg_ts = msg_ts + TD(0,0);
	// cout << msg_ts << endl;

	// Init para of matching
	arma::Col<double> stereo_TD_x(TD.size());
	arma::Col<double> stereo_TD_y(TD.size());
	arma::Col<double> stereo_TD_p(TD.size());
	arma::Col<double> stereo_TD_ts(TD.size());
	arma::Col<double> stereo_TD_rx(TD.size());
	arma::Col<double> stereo_TD_ry(TD.size());
	stereo_TD_x.fill(datum::nan);
	stereo_TD_y.fill(datum::nan);
	stereo_TD_p.fill(datum::nan);
	stereo_TD_ts.fill(datum::nan);
	stereo_TD_rx.fill(datum::nan);
	stereo_TD_ry.fill(datum::nan);

	arma::Col<double> stereo_TD_disparity_gt(TD.size());
	stereo_TD_disparity_gt.fill(datum::nan);

	// for each event
	for (int event_index = 1; event_index <= TD.n_rows; event_index++) {
		int current_event_index = event_index;
		double current_ts = TD(current_event_index - 1, TS_COL);
		int current_x = TD(current_event_index - 1, X_COL);
		int current_y = TD(current_event_index - 1, Y_COL);
		double current_p = TD(current_event_index - 1, P_COL);
		double current_disparity = TD(current_event_index - 1, GT_COL);
		if (TD(current_event_index - 1, C_COL) == LEFT_CAM) {
			arma::Col<double> min_error_row;
			int msg_data_bool = 0;
			for (int l_i = 0; l_i < LABEL; l_i++) {
				if (current_x - l_i > 1 && current_y - wradius > 1 && current_y + wradius < height) {
					arma::Col<double> sub_range = current_ts - last_spike_time_right(span(current_y - wradius - 1, current_y + wradius - 1), span(current_x - l_i - 1), span(0));
					arma::Col<double> sub_range_p = current_p - last_spike_time_right(span(current_y - wradius - 1, current_y + wradius - 1), span(current_x - l_i - 1), span(1));

					// consider the polarity
					uvec indexes = find(sub_range_p);
					Col<double> nan_vec(indexes.size());
					nan_vec.fill(datum::nan);
					sub_range(indexes) = nan_vec;
					
					// mark old events as invalid
					arma::Col<double> abs_sub_range = sub_range;
					abs_sub_range.transform([](double val) {return (val < 0 ? -1*val: val);});

					indexes = find(abs_sub_range > old_pixel_threshold);
					nan_vec.set_size(indexes.size());
					nan_vec.fill(datum::nan);
					sub_range(indexes) = nan_vec;

					// calculate the time difference
					abs_sub_range = sub_range;
					abs_sub_range.transform([](double val) {return (val < 0 ? -1*val: val);});
					Col<double> error_lr_t = abs_sub_range / temporal_var_Thre;

					// calculate distance to the epipolar line
					Col<double> error_lr_d = linspace(-wradius, wradius, wradius*2 + 1);
					error_lr_d.transform([](double val) {return (val < 0 ? -1*val: val);});
					error_lr_d = error_lr_d / special_var_Thre;
					Col<double> error_lr_row = error_lr_t + error_lr_d;
					indexes = find_finite(error_lr_row);

					// check if non nans have been found
					min_error_row.reset();
					if (!indexes.is_empty()) {
						min_error_row.set_size(1);
						min_error_row(0) = error_lr_row(indexes).min();
					}
				}

				if (!min_error_row.is_empty()) {
					msg_data(current_y - 1, current_x - 1, l_i) = min_error_row(0);
					msg_data_bool = msg_data_bool + 1;
				} else {
					msg_data(current_y - 1, current_x - 1, l_i) = DATE_TRUNC;
				}
			}
			if (current_ts == 153515) {
				int temp2 = 0;
			}
			if (msg_data_bool > 1) {
				if (current_y - WRADIUS_TS > 1 && current_y + WRADIUS_TS < height && current_x - WRADIUS_TS > 1 && current_x + WRADIUS_TS < width) {
					for (int i = 0; i < LABEL; i++) {
						Mat<double> bool_mrf_update = current_ts - msg_ts(span(current_y - WRADIUS_TS - 1, current_y + WRADIUS_TS - 1), span(current_x - WRADIUS_TS - 1, current_x + WRADIUS_TS - 1), span(i));
						bool_mrf_update.transform([](double val) {return (val < 0 ? -1*val: val);});
						bool_mrf_update.transform([](double val) {return (val < OLD_MRF_THRESHOLD);});

						Mat<double> msg_temp = msg_up(span(current_y - WRADIUS_TS - 1, current_y + WRADIUS_TS - 1), span(current_x - WRADIUS_TS - 1, current_x + WRADIUS_TS - 1), span(i));
						msg_up(span(current_y - WRADIUS_TS - 1, current_y + WRADIUS_TS - 1), span(current_x - WRADIUS_TS - 1, current_x + WRADIUS_TS - 1), span(i)) = msg_temp % bool_mrf_update;
						msg_temp = msg_down(span(current_y - WRADIUS_TS - 1, current_y + WRADIUS_TS - 1), span(current_x - WRADIUS_TS - 1, current_x + WRADIUS_TS - 1), span(i));
						msg_down(span(current_y - WRADIUS_TS - 1, current_y + WRADIUS_TS - 1), span(current_x - WRADIUS_TS - 1, current_x + WRADIUS_TS - 1), span(i)) = msg_temp % bool_mrf_update;
						msg_temp = msg_left(span(current_y - WRADIUS_TS - 1, current_y + WRADIUS_TS - 1), span(current_x - WRADIUS_TS - 1, current_x + WRADIUS_TS - 1), span(i));
						msg_left(span(current_y - WRADIUS_TS - 1, current_y + WRADIUS_TS - 1), span(current_x - WRADIUS_TS - 1, current_x + WRADIUS_TS - 1), span(i)) = msg_temp % bool_mrf_update;
						msg_temp = msg_right(span(current_y - WRADIUS_TS - 1, current_y + WRADIUS_TS - 1), span(current_x - WRADIUS_TS - 1, current_x + WRADIUS_TS - 1), span(i));
						msg_right(span(current_y - WRADIUS_TS - 1, current_y + WRADIUS_TS - 1), span(current_x - WRADIUS_TS - 1, current_x + WRADIUS_TS - 1), span(i)) = msg_temp % bool_mrf_update;
						msg_temp = msg_data(span(current_y - WRADIUS_TS - 1, current_y + WRADIUS_TS - 1), span(current_x - WRADIUS_TS - 1, current_x + WRADIUS_TS - 1), span(i));
						msg_data(span(current_y - WRADIUS_TS - 1, current_y + WRADIUS_TS - 1), span(current_x - WRADIUS_TS - 1, current_x + WRADIUS_TS - 1), span(i)) = msg_temp % bool_mrf_update;

					}
				}

				// bp iteration
				for (int iter = 1; iter <= BP_ITERATIONS; iter++) {
					// compute message
					// RIGHT
					for (int r = 0; r <= WRADIUS_TS; r++) {
						for (int r_y = -r; r_y <= r; r_y++) {
							for (int r_x = -r; r_x <= r; r_x++) {

								int y = current_y + r_y;
								int x = current_x + r_x;
								if (x > 1 && x < width && y > 1 && y < height) {
									// aggregate and find min
									double min_val = datum::inf;
									for (int value = 1; value <= LABEL; value++) {
										msg_left(y - 1, x, value - 1) = LAMBDA * msg_data(y - 1, x - 1, value - 1) + msg_left(y - 1, x - 1, value - 1) + msg_up(y - 1, x - 1, value - 1) + msg_down(y - 1, x - 1, value - 1);
										if (msg_left(y - 1, x, value - 1) < min_val) {
											min_val = msg_left(y - 1, x, value - 1);
										}
									}

									//dt
									double prev;
									for (int q = 2; q <= LABEL - 1; q++) {
										
										prev = msg_left(y - 1, x, q - 2) + 1.0;
										if (prev < msg_left(y - 1, x, q - 1)) {
											msg_left(y - 1, x, q - 1) = prev;
										}
									}

									for (int q = LABEL - 1; q >= 1; q--) {
										prev = msg_left(y - 1, x, q) + 1.0;
										if (prev < msg_left(y - 1, x, q - 1)) {
											msg_left(y - 1, x, q - 1) = prev;
										}
									}



									// truncate
									min_val = min_val + DISC_K; // tune DISC_K for smooth
									for (int value = 1; value <= LABEL; value++) {
										if (min_val < msg_left(y - 1, x, value - 1)) {
											msg_left(y - 1, x, value - 1) = min_val;
										}
									}

									// normalise
									double val = 0;
									for (int value = 1; value <= LABEL; value++) {
										val = val + msg_left(y - 1, x, value - 1);
									}

									val = val / LABEL;
									for (int value = 1; value <= LABEL; value++) {
										msg_left(y - 1, x, value - 1) = msg_left(y - 1, x, value - 1) - val;
									}


								}

							}
						}
					}
					// LEFT

					for (int r = 0; r <= WRADIUS_TS; r++) {
						for (int r_y = -r; r_y <= r; r_y++) {
							for (int r_x = -r; r_x <= r; r_x++) {
								int y = current_y + r_y;
								int x = current_x + r_x;
								if (x > 1 && x < width && y > 1 && y < height) {
									// aggregate and find min
									double min_val = datum::inf;
									for (int value = 1; value <= LABEL; value++) {
										msg_right(y - 1, x - 2, value - 1) = LAMBDA * msg_data(y - 1, x - 1, value - 1) + msg_right(y - 1, x - 1, value - 1) + msg_up(y - 1, x - 1, value - 1) + msg_down(y - 1, x - 1, value - 1);
										if (msg_right(y - 1, x - 2, value - 1) < min_val) {
											min_val = msg_right(y - 1, x - 2, value - 1);
										}
									}
									//dt
									double prev;
									for (int q = 2; q <= LABEL - 1; q++) {
										prev = msg_right(y - 1, x - 2, q - 2) + 1.0;
										if (prev < msg_right(y - 1, x - 2, q - 1)) {
											msg_right(y - 1, x - 2, q - 1) = prev;
										}
									}
									for (int q = LABEL - 1; q >= 1; q--) {
										prev = msg_right(y - 1, x - 2, q) + 1.0;
										if (prev < msg_right(y - 1, x - 2, q - 1)) {
											msg_right(y - 1, x - 2, q - 1) = prev;
										}
									}
									// truncate
									min_val = min_val + DISC_K; // tune DISC_K for smooth
									for (int value = 1; value <= LABEL; value++) {
										if (min_val < msg_right(y - 1, x - 2, value - 1)) {
											msg_right(y - 1, x - 2, value - 1) = min_val;
										}
									}
									// normalise
									double val = 0;
									for (int value = 1; value <= LABEL; value++) {
										val = val + msg_right(y - 1, x - 2, value - 1);
									}
									val = val / LABEL;
									for (int value = 1; value <= LABEL; value++) {
										msg_right(y - 1, x - 2, value - 1) = msg_right(y - 1, x - 2, value - 1) - val;
									}
								}
							}
						}
					}
					// DOWN
					for (int r = 0; r <= WRADIUS_TS; r++) {
						for (int r_y = -r; r_y <= r; r_y++) {
							for (int r_x = -r; r_x <= r; r_x++) {
								int y = current_y + r_y;
								int x = current_x + r_x;
								if (x > 1 && x < width && y > 1 && y < height) {
									// aggregate and find min
									double min_val = datum::inf;
									for (int value = 1; value <= LABEL; value++) {
										msg_up(y, x - 1, value - 1) = LAMBDA * msg_data(y - 1, x - 1, value - 1) + msg_left(y - 1, x - 1, value - 1) + msg_up(y - 1, x - 1, value - 1) + msg_right(y - 1, x - 1, value - 1);
										if (msg_up(y, x - 1, value - 1) < min_val) {
											min_val = msg_up(y, x - 1, value - 1);
										}
									}
									//dt
									double prev;
									for (int q = 2; q <= LABEL - 1; q++) {
										prev = msg_up(y, x - 1, q - 2) + 1.0;
										if (prev < msg_up(y, x - 1, q - 1)) {
											msg_up(y, x - 1, q - 1) = prev;
										}
									}
									for (int q = LABEL - 1; q >= 1; q--) {
										prev = msg_up(y, x - 1, q) + 1.0;
										if (prev < msg_up(y, x - 1, q - 1)) {
											msg_up(y, x - 1, q - 1) = prev;
										}
									}
									// truncate
									min_val = min_val + DISC_K; // tune DISC_K for smooth
									for (int value = 1; value <= LABEL; value++) {
										if (min_val < msg_up(y, x - 1, value - 1)) {
											msg_up(y, x - 1, value - 1) = min_val;
										}
									}
									// normalise
									double val = 0;
									for (int value = 1; value <= LABEL; value++) {
										val = val + msg_up(y, x - 1, value - 1);
									}
									val = val / LABEL;
									for (int value = 1; value <= LABEL; value++) {
										msg_up(y, x - 1, value - 1) = msg_up(y, x - 1, value - 1) - val;
									}
								}
							}
						}
					}
					// UP
					for (int r = 0; r <= WRADIUS_TS; r++) {
						for (int r_y = -r; r_y <= r; r_y++) {
							for (int r_x = -r; r_x <= r; r_x++) {
								int y = current_y + r_y;
								int x = current_x + r_x;
								if (x > 1 && x < width && y > 1 && y < height) {
									// aggregate and find min
									double min_val = datum::inf;
									for (int value = 1; value <= LABEL; value++) {
										msg_down(y - 2, x - 1, value - 1) = LAMBDA * msg_data(y - 1, x - 1, value - 1) + msg_left(y - 1, x - 1, value - 1) + msg_right(y - 1, x - 1, value - 1) + msg_down(y - 1, x - 1, value - 1);
										if (msg_down(y - 2, x - 1, value - 1) < min_val) {
											min_val = msg_down(y - 2, x - 1, value - 1);
										}
									}
									//dt
									double prev;
									for (int q = 2; q <= LABEL - 1; q++) {
										prev = msg_down(y - 2, x - 1, q - 2) + 1.0;
										if (prev < msg_down(y - 2, x - 1, q - 1)) {
											msg_down(y - 2, x - 1, q - 1) = prev;
										}
									}
									for (int q = LABEL - 1; q >= 1; q--) {
										prev = msg_down(y - 2, x - 1, q) + 1.0;
										if (prev < msg_down(y - 2, x - 1, q - 1)) {
											msg_down(y, x - 1, q - 1) = prev;
										}
									}
									// truncate
									min_val = min_val + DISC_K; // tune DISC_K for smooth
									for (int value = 1; value <= LABEL; value++) {
										if (min_val < msg_down(y - 2, x - 1, value - 1)) {
											msg_down(y - 2, x - 1, value - 1) = min_val;
										}
									}
									// normalise
									double val = 0;
									for (int value = 1; value <= LABEL; value++) {
										val = val + msg_down(y - 2, x - 1, value - 1);
									}
									val = val / LABEL;
									for (int value = 1; value <= LABEL; value++) {
										msg_down(y - 2, x - 1, value - 1) = msg_down(y - 2, x - 1, value - 1) - val;
									}
								}
							}
						}
					}
				}
				// updatetime
				if (current_y - WRADIUS_TS >= 1 && current_y + WRADIUS_TS <= height && current_x - WRADIUS_TS >= 1 && current_x + WRADIUS_TS <= width) {
					msg_ts(span(current_y - WRADIUS_TS - 1, current_y + WRADIUS_TS - 1), span(current_x - WRADIUS_TS - 1, current_x + WRADIUS_TS - 1), span::all).fill(current_ts);
				}
				// Finds the MAP assignment as well as calculating energy
				// MAP assignment
				int y = current_y;
				int x = current_x;
				double best = datum::inf;
				for (int j = 1; j <= LABEL; j++) {
					double cost_event = msg_left(y - 1, x - 1, j - 1) + msg_right(y - 1, x - 1, j - 1) + msg_up(y - 1, x - 1, j - 1) + msg_down(y - 1, x - 1, j - 1) + LAMBDA * msg_data(y - 1, x - 1, j - 1);
					if (cost_event < best) {
						best = cost_event;
						disparity = j - 1;
					}
				}
				if (best < THETA && disparity != 0) {
					stereo_TD_x(event_index - 1) = current_x;
					stereo_TD_y(event_index - 1) = current_y;
					stereo_TD_p(event_index - 1) = disparity;
					stereo_TD_ts(event_index - 1) = current_ts;
					if (groundtruth_bool) {
						stereo_TD_disparity_gt(event_index - 1) = current_disparity;
					}
				}
			}

		} else {
			last_spike_time_right(TD(event_index - 1, Y_COL) - 1, TD(event_index - 1, X_COL) - 1, 0) = TD(event_index - 1, TS_COL);
			last_spike_time_right(TD(event_index - 1, Y_COL) - 1, TD(event_index - 1, X_COL) - 1, 1) = TD(event_index - 1, P_COL);
		}

	}

	uvec indexes = find_finite(stereo_TD_ts);
	stereo_TD_x = stereo_TD_x(indexes);
	stereo_TD_y = stereo_TD_y(indexes);
	stereo_TD_p = stereo_TD_p(indexes);
	stereo_TD_ts = stereo_TD_ts(indexes);

	if (groundtruth_bool) {
		stereo_TD_disparity_gt = stereo_TD_disparity_gt(indexes);
	}

	Mat<double> stereo_TD;
	if (groundtruth_bool) {
		stereo_TD.set_size(stereo_TD_x.size(), 5);
		stereo_TD(span::all, span(4)) = stereo_TD_disparity_gt;	
	} else {
		stereo_TD.set_size(stereo_TD_x.size(), 4);
	}
	stereo_TD(span::all, span(0)) = stereo_TD_x;
	stereo_TD(span::all, span(1)) = stereo_TD_y;
	stereo_TD(span::all, span(2)) = stereo_TD_p;
	stereo_TD(span::all, span(3)) = stereo_TD_ts;

	return stereo_TD;
}
