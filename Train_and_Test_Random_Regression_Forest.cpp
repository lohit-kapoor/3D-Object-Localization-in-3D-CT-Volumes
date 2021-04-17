#include <iostream>
#include<conio.h>
#include<math.h>
#include<process.h>
#include <fstream>
#include <string>
#include <sstream>
#include <stdlib.h>
#include <omp.h>
#include <vector>
#include<algorithm>
#include <queue>
#include <map>
#include <Windows.h>

typedef std::vector<std::vector<long double> > vvld;
#define totrows 845864


class Node
{
public:
	int split_feature_id;
	long int split_feature_value;
	long double tr_cov;
	long double d_mean[6];
	int depth_node;
	bool isLeaf;
	vvld leaf_test_data;
	Node* left;
	Node* right;


};

std::map<long double, Node> forest_leaves;

long double det(long double b[6][6], int m)
{
	int i, j;
	long double sum = 0.0;
	long double c[6][6] = { 0 };
	if (m == 2)
	{
		sum = b[0][0] * b[1][1] - b[0][1] * b[1][0];
		return sum;
	}
	for (int p = 0; p<m; p++)
	{
		int h = 0, k = 0;
		for (i = 1; i<m; i++)
		{
			for (j = 0; j<m; j++)
			{
				if (j == p)
				{
					continue;
				}
				c[h][k] = b[i][j];
				k++;
				if (k == m - 1)
				{
					h++;
					k = 0;
				}

			}
		}

		sum = sum + b[0][p] * pow(-1, p)*det(c, m - 1);
	}
	return sum;
}



long double find_det_of_cov(vvld& fdata)
{
	long double sum_d1 = 0, sum_d2 = 0, sum_d3 = 0, sum_d4 = 0, sum_d5 = 0, sum_d6 = 0;
	for (unsigned int p = 0; p < fdata.size(); p++)
	{
		sum_d1 += fdata[p][0];
		sum_d2 += fdata[p][1];
		sum_d3 += fdata[p][2];
		sum_d4 += fdata[p][3];
		sum_d5 += fdata[p][4];
		sum_d6 += fdata[p][5];
	}

	long double mean_d1 = (long double)sum_d1 / (long double)fdata.size();
	long double mean_d2 = (long double)sum_d2 / (long double)fdata.size();
	long double mean_d3 = (long double)sum_d3 / (long double)fdata.size();
	long double mean_d4 = (long double)sum_d4 / (long double)fdata.size();
	long double mean_d5 = (long double)sum_d5 / (long double)fdata.size();
	long double mean_d6 = (long double)sum_d6 / (long double)fdata.size();
	long double mean[6] = { mean_d1, mean_d2, mean_d3, mean_d4, mean_d5, mean_d6 };

	vvld Y(fdata.size(), std::vector<long double>(6));
	vvld YT(6, std::vector<long double>(fdata.size()));

	for (unsigned int g = 0; g < fdata.size(); g++)
	{
		for (unsigned int h = 0; h < 6; h++)
		{
			Y[g][h] = fdata[g][h] - mean[h];
			YT[h][g] = Y[g][h];
		}
	}

	long double cov[6][6] = { 0 };
	for (int i = 0; i < 6; ++i)
	{
		for (int j = 0; j < 6; ++j)
		{
			for (unsigned int k = 0; k < fdata.size(); ++k)
			{
				cov[i][j] += YT[i][k] * Y[k][j];

			}
			(long double)cov[i][j] = (long double)cov[i][j] / ((long double)fdata.size() - 1);
		}
	}
	Y.clear(); YT.clear();

	return det(cov, 6);

}



Node* buildDecisionTree(Node* Node_ptr, vvld& data, int Node_depth, long double det_data)
{
	double end3 = omp_get_wtime();
	Node_ptr->depth_node = Node_depth;
	std::cout << "Det_Cov_NodeData: " << det_data << std::endl;

	if (Node_depth < 8 && data.size() >100)
	{
		Node_ptr->isLeaf = false;
		int feature_id = rand() % 10;
		std::cout << "feature_id: " << feature_id << std::endl;
		Node_ptr->split_feature_id = feature_id;


		long double feature_min = data[0][feature_id + 6];
		long double feature_max = 0;
		for (unsigned int p = 0; p < data.size(); p++)
		{
			if (data[p][feature_id + 6] > feature_max)
				feature_max = data[p][feature_id + 6];
			if (data[p][feature_id + 6] < feature_min)
				feature_min = data[p][feature_id + 6];
		}


		long int range = (long int)feature_max - (long int)feature_min;
		std::cout << "range: " << range << " feature_max: " << feature_max << " feature_min: " << feature_min << std::endl;
		std::vector<long int> fsplit(10, 0);
		std::vector<long double> IG(10, 0);  int ii;
		std::vector<long double> det_cov_rightdata(10, 0);
		std::vector<long double> det_cov_leftdata(10, 0);

		for (ii = 0; ii < 10; ii++)
		{
			fsplit[ii] = rand() % range + (long int)feature_min;
			std::cout << "fsplit[" << ii << "]: " << fsplit[ii] << std::endl;
			vvld temp_left_data;
			vvld temp_right_data;
			for (unsigned int gg = 0; gg < data.size(); gg++)
			{
				std::vector<long double> row(16);
				for (int dd = 0; dd < 16; dd++)
				{
					row[dd] = data[gg][dd];
				}
				if (data[gg][feature_id + 6] >= fsplit[ii])
				{
					temp_right_data.push_back(row);
				}
				else
				{
					temp_left_data.push_back(row);
				}
				row.clear();
			}
			if (temp_right_data.size() == 1 || temp_right_data.size() == 0){ det_cov_rightdata[ii] = 0; }
			else
			{
				det_cov_rightdata[ii] = find_det_of_cov(temp_right_data);
			}
			if (temp_left_data.size() == 1 || temp_left_data.size() == 0){ det_cov_leftdata[ii] = 0; }
			else
			{
				det_cov_leftdata[ii] = find_det_of_cov(temp_left_data);
			}
			std::cout << " data_size: " << data.size() << " temp_left_datasize: " << temp_left_data.size() << " temp_right_datasize: " << temp_right_data.size() << std::endl;
			long double w_left = (long double)temp_left_data.size() / (long double)data.size();
			long double w_right = (long double)temp_right_data.size() / (long double)data.size();
			std::cout << " det_cov_rightdata[" << ii << "]: " << det_cov_rightdata[ii] << " det_cov_leftdata[" << ii << "]: " << det_cov_leftdata[ii] << " w_left: " << w_left << "w_right: " << w_right << std::endl;
			IG[ii] = log(1 + std::abs(det_data)) - (w_left * log(1 + std::abs(det_cov_leftdata[ii]))) - (w_right * log(1 + std::abs(det_cov_rightdata[ii])));
			std::cout << "IG[" << ii << "]: " << IG[ii] << std::endl;
			temp_left_data.clear(); temp_right_data.clear();
		}

		auto p_max_IG = std::max_element(std::begin(IG), std::end(IG));
		long double max_IG = *p_max_IG;
		std::cout << "max_IG: " << max_IG << std::endl;
		if (max_IG <= 0) goto shifthere;
		int max_index = std::distance(std::begin(IG), p_max_IG);
		long int max_fsplit = fsplit[max_index];
		long double max_det_right = det_cov_rightdata[max_index];
		long double max_det_left = det_cov_leftdata[max_index];
		Node_ptr->split_feature_value = max_fsplit;
		fsplit.clear(); IG.clear(); det_cov_rightdata.clear(); det_cov_leftdata.clear();

		vvld left_data, right_data;
		for (unsigned int hh = 0; hh < data.size(); hh++)
		{
			std::vector<long double> row(16);
			for (int dd = 0; dd < 16; dd++)
			{
				row[dd] = data[hh][dd];
			}
			if (data[hh][feature_id + 6] >= max_fsplit)
			{
				right_data.push_back(row);
			}
			else
			{
				left_data.push_back(row);
			}
			row.clear();

		}


		Node* right_child = new Node;
		Node* left_child = new Node;
		double end4 = omp_get_wtime();
		std::cout << "Done building node at depth: " << Node_depth << " with IG: " << max_IG << " in time: " << end4 - end3 << " seconds" << std::endl;

		Node_depth++; data.clear();
		Node_ptr->right = buildDecisionTree(right_child, right_data, Node_depth, max_det_right);
		Node_ptr->left = buildDecisionTree(left_child, left_data, Node_depth, max_det_left);


	}

	else
	{
	shifthere:
		Node_ptr->isLeaf = true;
		long double sum_d1 = 0, sum_d2 = 0, sum_d3 = 0, sum_d4 = 0, sum_d5 = 0, sum_d6 = 0;
		for (unsigned int p = 0; p < data.size(); p++)
		{
			sum_d1 += data[p][0];
			sum_d2 += data[p][1];
			sum_d3 += data[p][2];
			sum_d4 += data[p][3];
			sum_d5 += data[p][4];
			sum_d6 += data[p][5];
		}

		long double mean_d1 = (long double)sum_d1 / (long double)data.size();
		long double mean_d2 = (long double)sum_d2 / (long double)data.size();
		long double mean_d3 = (long double)sum_d3 / (long double)data.size();
		long double mean_d4 = (long double)sum_d4 / (long double)data.size();
		long double mean_d5 = (long double)sum_d5 / (long double)data.size();
		long double mean_d6 = (long double)sum_d6 / (long double)data.size();
		long double mean[6] = { mean_d1, mean_d2, mean_d3, mean_d4, mean_d5, mean_d6 };
		for (int u = 0; u < 6; u++)
			Node_ptr->d_mean[u] = mean[u];

		vvld Y(data.size(), std::vector<long double>(6));
		vvld YT(6, std::vector<long double>(data.size()));

		for (unsigned int g = 0; g < data.size(); g++)
		{
			for (unsigned int h = 0; h < 6; h++)
			{
				Y[g][h] = data[g][h] - mean[h];
				YT[h][g] = Y[g][h];
			}
		}

		vvld cov(6, std::vector<long double>(6));
		for (int i = 0; i < 6; ++i)
		{
			for (int j = 0; j < 6; ++j)
			{
				for (unsigned int k = 0; k < data.size(); ++k)
				{
					cov[i][j] += YT[i][k] * Y[k][j];

				}
				(long double)cov[i][j] = (long double)cov[i][j] / ((long double)data.size());
			}
		}
		Y.clear(); YT.clear();
		long double trace = 0;

		for (int j = 0; j < 6; ++j)
		{
			trace += cov[j][j];
		}

		Node_ptr->tr_cov = trace;
		Node_ptr->right = NULL;
		Node_ptr->left = NULL;
		cov.clear(); data.clear();
		double end4 = omp_get_wtime();
		std::cout << "Done building leaf at depth: " << Node_depth << " in time: " << end4 - end3 << " seconds with cov_trace: " << trace << std::endl;
		return Node_ptr;


	}


	return Node_ptr;
}

void testForest(Node* ptr_Node, vvld& test_data)
{
	if (ptr_Node->isLeaf)
	{
		if (ptr_Node->tr_cov != 0)
		{
			ptr_Node->leaf_test_data = test_data;
			forest_leaves[ptr_Node->tr_cov] = *ptr_Node;
			std::cout << "Test Data reaching this leaf is stored in it with trace: "<< ptr_Node->tr_cov << " and test data stored size: "<< test_data.size() <<std::endl;
		}
		test_data.clear();
	}
	else
	{
		vvld right_data, left_data;
		for (unsigned int hh = 0; hh < test_data.size(); hh++)
		{
			std::vector<long double> row(16);
			for (int dd = 0; dd < 16; dd++)
			{
				row[dd] = test_data[hh][dd];
			}
			if (test_data[hh][ptr_Node->split_feature_id + 6] >= ptr_Node->split_feature_value)
			{
				right_data.push_back(row);
			}
			else
			{
				left_data.push_back(row);
			}
			row.clear();

		}
		test_data.clear();
		if (!right_data.empty()) testForest(ptr_Node->right, right_data);
		if (!left_data.empty()) testForest(ptr_Node->left, left_data);


	}

}


void Levelorder(Node* root_node)
{
	std::queue<Node*> Q;
	Q.push(root_node);
	while (!Q.empty())
	{
		Node* Current = Q.front();
		std::cout << " Depth: " << Current->depth_node << "Node is leaf?: " << Current->isLeaf << std::endl;
		if (Current->right != NULL) Q.push(Current->right);
		if (Current->left != NULL) Q.push(Current->left);
		Q.pop();

	}

}

void LevelOrderLeafDataRemoval(Node* root_node)
{
	std::queue<Node*> Q2;
	Q2.push(root_node);
	while (!Q2.empty())
	{
		Node* Current = Q2.front();
		if (Current->isLeaf)
		{
			if (!Current->leaf_test_data.empty()) Current->leaf_test_data.clear();
		}
		if (Current->right != NULL) Q2.push(Current->right);
		if (Current->left != NULL) Q2.push(Current->left);
		Q2.pop();
	}
}

int main()
{
	Node* root_0 = new Node;
	Node* root_1 = new Node;
	Node* root_2 = new Node;
	Node* root_3 = new Node;
	Node* root_4 = new Node;
	Node* root_5 = new Node;
	Node* root_6 = new Node;
	Node* root_7 = new Node;
	Node* root_8 = new Node;
	Node* root_9 = new Node;
	Node* root_10 = new Node;
	Node* root_11 = new Node;

	double start = omp_get_wtime();

	vvld train_data(totrows, std::vector<long double>(16));

	std::ifstream file("E:\\Train_Data_226.csv");
	std::string line;
	int col = 0;
	int row = 0;
	while (std::getline(file, line))
	{
		std::istringstream iss(line);
		std::string result;
		while (std::getline(iss, result, ','))
		{
			train_data[row][col] = std::stold(result.c_str());
			col = col + 1;
		}
		row = row + 1;
		col = 0;
	}

	file.close();
	double end1 = omp_get_wtime();
	std::cout << "Done Uploading Training Data in time: " << end1 - start << "seconds" << std::endl;
	long double det_cov_traindata = find_det_of_cov(train_data);

	double end2 = omp_get_wtime();
	std::cout << "Done finding det of cov of Training Data in time: " << end2 - end1 << "seconds" << std::endl;

	for (int trees = 0; trees < 12; trees++)
	{

		vvld tree_train_data = train_data;

		int depth = 0;

		double end22 = omp_get_wtime();

		if (trees == 0) root_0 = buildDecisionTree(root_0, tree_train_data, depth, det_cov_traindata);
		else if (trees == 1) root_1 = buildDecisionTree(root_1, tree_train_data, depth, det_cov_traindata);
		else if (trees == 2) root_2 = buildDecisionTree(root_2, tree_train_data, depth, det_cov_traindata);
		else if (trees == 3) root_3 = buildDecisionTree(root_3, tree_train_data, depth, det_cov_traindata);
		else if (trees == 4) root_4 = buildDecisionTree(root_4, tree_train_data, depth, det_cov_traindata);
		else if (trees == 5) root_5 = buildDecisionTree(root_5, tree_train_data, depth, det_cov_traindata);
		else if (trees == 6) root_6 = buildDecisionTree(root_6, tree_train_data, depth, det_cov_traindata);
		else if (trees == 7) root_7 = buildDecisionTree(root_7, tree_train_data, depth, det_cov_traindata);
		else if (trees == 8) root_8 = buildDecisionTree(root_8, tree_train_data, depth, det_cov_traindata);
		else if (trees == 9) root_9 = buildDecisionTree(root_9, tree_train_data, depth, det_cov_traindata);
		else if (trees == 10) root_10 = buildDecisionTree(root_10, tree_train_data, depth, det_cov_traindata);
		else if (trees == 11) root_11 = buildDecisionTree(root_11, tree_train_data, depth, det_cov_traindata);

		double end112 = omp_get_wtime();
		std::cout << "Time taken to build one Decision Tree: " << end112 - end22 << " seconds" << std::endl;

	}
	train_data.clear();


	Levelorder(root_0);
	Levelorder(root_1);
	Levelorder(root_2);
	Levelorder(root_3);
	Levelorder(root_4);
	Levelorder(root_5);
	Levelorder(root_6);
	Levelorder(root_7);
	Levelorder(root_8);
	Levelorder(root_9);
	Levelorder(root_10);
	Levelorder(root_11);

	
	for (int testing = 0; testing < 30; testing++)
	{
		Beep(523, 5000);
		std::string s_datapoints, s_filename, s_ntrees;
		unsigned int datapoints; int ntrees;
		std::cin >> s_datapoints;
		std::cin >> s_filename;
		std::cin >> s_ntrees;
		datapoints = std::stoul(s_datapoints.c_str());
		ntrees = std::stoi(s_ntrees.c_str());
		unsigned int percentdata = (unsigned int)datapoints / (unsigned int)100;
		std::cout <<"Number of DataPoints entered: " << datapoints << std::endl;
		std::cout << "Filename entered" << s_filename << std::endl;
		std::cout << "Number of Trees entered" << ntrees << std::endl;
		
		vvld test_data(datapoints, std::vector<long double>(16));
		std::ifstream test_file(s_filename);
		std::string test_line;
		int test_col = 0;
		int test_row = 0;
		while (std::getline(test_file, test_line))
		{
			std::istringstream iss(test_line);
			std::string test_result;
			while (std::getline(iss, test_result, ','))
			{
				test_data[test_row][test_col] = std::stold(test_result.c_str());
				test_col = test_col + 1;
			}
			test_row = test_row + 1;
			test_col = 0;
		}

		test_file.close();




		for (int trees = 0; trees < ntrees; trees++)
		{

			vvld tree_test_data = test_data;

			double end_t1_1 = omp_get_wtime();

			if (trees == 0) testForest(root_0, tree_test_data);
			else if (trees == 1) testForest(root_1, tree_test_data);
			else if (trees == 2) testForest(root_2, tree_test_data);
			else if (trees == 3) testForest(root_3, tree_test_data);
			else if (trees == 4) testForest(root_4, tree_test_data);
			else if (trees == 5) testForest(root_5, tree_test_data);
			else if (trees == 6) testForest(root_6, tree_test_data);
			else if (trees == 7) testForest(root_7, tree_test_data);
			else if (trees == 8) testForest(root_8, tree_test_data);
			else if (trees == 9) testForest(root_9, tree_test_data);
			else if (trees == 10) testForest(root_10, tree_test_data);
			else if (trees == 11) testForest(root_11, tree_test_data);

			double end_t1_2 = omp_get_wtime();
			std::cout << "Time taken to test one Decision Tree: " << end_t1_2 - end_t1_1 << " seconds" << std::endl;

		}

		test_data.clear();
		vvld bounding_box;
		int voxelNumber = 0;
		for (std::map<long double, Node>::iterator it = forest_leaves.begin(); it != forest_leaves.end(); it++)
		{
			std::cout << "This testing leaf has Trace: " << it->first << " with Data Size: " << it->second.leaf_test_data.size() << std::endl;

			for (unsigned int pp = 0; pp < it->second.leaf_test_data.size(); pp++)
			{

				voxelNumber++;
				std::cout << "Voting Voxel Number: " << voxelNumber << std::endl;
				std::vector<long double> bbrow(6);
				for (int f = 0; f < 6; f++)
				{
					bbrow[f] = it->second.leaf_test_data[pp][f] - it->second.d_mean[f];
					std::cout << " has Location in mm: " << it->second.leaf_test_data[pp][f] << ", ";
				}std::cout << std::endl;
				if (bounding_box.size() < percentdata)
				{
					bounding_box.push_back(bbrow);
					bbrow.clear();
				}
				else goto exithere;
			}
			
		}

	exithere:
		long double sum_b1 = 0, sum_b2 = 0, sum_b3 = 0, sum_b4 = 0, sum_b5 = 0, sum_b6 = 0;
		for (unsigned int p = 0; p < bounding_box.size(); p++)
		{
			sum_b1 += bounding_box[p][0];
			sum_b2 += bounding_box[p][1];
			sum_b3 += bounding_box[p][2];
			sum_b4 += bounding_box[p][3];
			sum_b5 += bounding_box[p][4];
			sum_b6 += bounding_box[p][5];
		}

		long double mean_b1 = (long double)sum_b1 / (long double)bounding_box.size();
		long double mean_b2= (long double)sum_b2 / (long double)bounding_box.size();
		long double mean_b3 = (long double)sum_b3 / (long double)bounding_box.size();
		long double mean_b4= (long double)sum_b4/ (long double)bounding_box.size();
		long double mean_b5 = (long double)sum_b5 / (long double)bounding_box.size();
		long double mean_b6 = (long double)sum_b6/ (long double)bounding_box.size();
		long double mean_bb[6] = { mean_b1, mean_b2, mean_b3, mean_b4, mean_b5, mean_b6 };
		std::cout << "Absolute Bounding Box Prediction for Test Data : " << std::endl;
		for (int u = 0; u < 6; u++)
			std::cout << "BB[" << u << "]" << mean_bb[u] << ",";
		std::cout << std::endl;
		forest_leaves.clear(); bounding_box.clear();

		for (int rtrees = 0; rtrees < ntrees; rtrees++)
		{
			if (rtrees == 0) LevelOrderLeafDataRemoval(root_0);
			else if (rtrees == 1) LevelOrderLeafDataRemoval(root_1);
			else if (rtrees == 2) LevelOrderLeafDataRemoval(root_2);
			else if (rtrees == 3) LevelOrderLeafDataRemoval(root_3);
			else if (rtrees == 4) LevelOrderLeafDataRemoval(root_4);
			else if (rtrees == 5) LevelOrderLeafDataRemoval(root_5);
			else if (rtrees == 6) LevelOrderLeafDataRemoval(root_6);
			else if (rtrees == 7) LevelOrderLeafDataRemoval(root_7);
			else if (rtrees == 8) LevelOrderLeafDataRemoval(root_8);
			else if (rtrees == 9) LevelOrderLeafDataRemoval(root_9);
			else if (rtrees == 10) LevelOrderLeafDataRemoval(root_10);
			else if (rtrees == 11) LevelOrderLeafDataRemoval(root_11);
			
		}


	}
	return 0;

}