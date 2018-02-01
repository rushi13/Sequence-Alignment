#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdlib.h> 
#include <vector>
#include <algorithm>
#include <ctype.h>

using namespace std;

struct alignment
{
	string query_id;
	string database_id;
	string query_sequence;
	string database_sequence;
	int query_start;
	int database_start;
	int query_end;
	int database_end;
	int score;
};

struct FirstColumnOnlyCmp
{
    bool operator()(const alignment& x, const alignment& y) const
    {
        return x.score > y.score; 
    }
};

std::string to_string(int i)
{
    std::stringstream ss;
    ss << i;
    return ss.str();
}

int direction()
{
	int dir = -1;
	
	return dir;
}

int score(string alphabet,vector <vector <int> > &score_mat, char query, char database)
{
	int i = 0;

	while(tolower(alphabet[i]) != tolower(query))
		i++;

	int j = 0;

	while(tolower(alphabet[j]) != tolower(database))
		j++;

	return score_mat[i][j];
}

alignment dovetail_alignment(int penalty,string query_id,string database_id, string query, string database,string alphabet,vector <vector <int> > &score_mat)
{
	int result = 0,query_start = -1, query_end = -1, database_start = 1, database_end = -1, dt = 0;
	
	string laq = "", lad = "",tempd = "",tempq = "";
	
	int **dp = 0;
    
	dp = new int *[query.length() + 1];
	
	for( int i = 0 ; i < query.length() + 1 ; i++ )
        dp[i] = new int [database.length() + 1];
	
	int **bt = 0;

    bt = new int *[query.length() + 1];
	
	for( int i = 0 ; i < query.length() + 1 ; i++ )
        bt[i] = new int [database.length() + 1];

	int i = 0, j = 0;
	
	for(i = 0;i <= query.length();i++)
		dp[i][0] = 0;
	
	for(i = 0;i <= database.length();i++)
		dp[0][i] = 0;
	
	for(i = 0;i <= query.length();i++)
		bt[i][0] = 2;

	for(i = 0;i <= database.length();i++)
		bt[0][i] = 3;
	
	for(int i = 1;i <= query.length(); i++)
	{
		for(int j = 1;j <= database.length();j++)
		{
			dp[i][j] = max(dp[i-1][j-1] + score(alphabet,score_mat,query[i-1],database[j-1]),max(dp[i-1][j] + penalty,dp[i][j-1] + penalty));
			
			if(max(dp[i-1][j-1] + score(alphabet,score_mat,query[i-1],database[j-1]),max(dp[i-1][j] + penalty,dp[i][j-1] + penalty)) == (dp[i-1][j-1] + score(alphabet,score_mat,query[i-1],database[j-1])))
				bt[i][j] = 1;
			else if(max(dp[i-1][j-1] + score(alphabet,score_mat,query[i-1],database[j-1]),max(dp[i-1][j] + penalty,dp[i][j-1] + penalty)) == (dp[i-1][j] + penalty))
				bt[i][j] = 2;
			else
				bt[i][j] = 3;
			
		}
	}
		
	for(j = 1;j < database.length() + 1;j++)
	{
		if(dp[query.length()][j] > result)
		{
			result = dp[query.length()][j];
			database_end = j;
			database_start = 0;
			query_end = query.length() - 1;
			dt = 1;
		}
	}


	for(i = 1;i < query.length() + 1;i++)
	{
		if(dp[i][database.length()] > result)
		{
			result = dp[i][database.length()];
			query_end = i;
			query_start = 0;
			database_start = -1;
			database_end = database.length() - 1;
			dt = 2;
		}
	}
	
			
	if(dt == 1)
	{
		j = database_end;
		i = query.length();
	}
	else
	{
		i = query_end;
		j = database.length();
	}
	
	while(i > 0 && j > 0)
	{
		tempd = "";
		tempq = "";
		if(bt[i][j] == 1)
		{
			i--;
			j--;
			tempq.push_back(query[i]);
			tempd.push_back(database[j]);
			laq.insert(0,tempq);
			lad.insert(0,tempd);
		}
		else if(bt[i][j] == 2)
		{
			i--;
			tempq.push_back(query[i]);
			laq.insert(0,tempq);
			lad.insert(0,".");
		}
		else if(bt[i][j] == 3)
		{
			j--;
			tempd.push_back(database[j]);
			laq.insert(0,".");
			lad.insert(0,tempd);
		}
		
	}
	
	if(dt == 1)
	{
		query_start = i;
	}
	else
	{
		database_start = j;
	}
	
	alignment res;
	
	res.query_id = query_id;
	res.database_id = database_id;
	res.query_sequence = laq;
	res.database_sequence = lad;
	res.query_start = query_start;
	res.database_start = database_start;
	res.query_end = query_end;
	res.database_end = database_end;
	res.score = result;
	
	for(int i = 0; i < query.length() + 1; ++i) {
        delete[] bt[i];   
    }
    
    delete[] bt;
	
	for(int i = 0; i < query.length() + 1; ++i) {
        delete[] dp[i];   
    }

    delete[] dp;
	
	return res;	
}

alignment local_alignment(int penalty,string query_id,string database_id,string query, string database,string alphabet,vector <vector <int> > &score_mat)
{
	int result = 0,query_start = -1, query_end = -1, database_start = -1, database_end = -1;
	
	string laq = "", lad = "",tempd = "",tempq = "";
	
	int **dp = 0;
    
	dp = new int *[query.length() + 1];
	
	for( int i = 0 ; i < query.length() + 1 ; i++ )
        dp[i] = new int [database.length() + 1];
	
	int **bt = 0;

    bt = new int *[query.length() + 1];
	
	for( int i = 0 ; i < query.length() + 1 ; i++ )
        bt[i] = new int [database.length() + 1];

	int i = 0, j = 0;
	
	for(i = 0;i <= query.length();i++)
		dp[i][0] = 0;
	
	for(i = 0;i <= database.length();i++)
		dp[0][i] = 0;
	
	for(i = 0;i <= query.length();i++)
		bt[i][0] = 2;

	for(i = 0;i <= database.length();i++)
		bt[0][i] = 3;
	
	for(int i = 1;i <= query.length(); i++)
	{
		for(int j = 1;j <= database.length();j++)
		{
			dp[i][j] = max(0,max(dp[i-1][j-1] + score(alphabet,score_mat,query[i-1],database[j-1]),max(dp[i-1][j] + penalty,dp[i][j-1] + penalty)));
			
			if(dp[i][j] == (dp[i-1][j-1] + score(alphabet,score_mat,query[i-1],database[j-1])))
				bt[i][j] = 1;
			else if(dp[i][j] == (dp[i-1][j] + penalty))
				bt[i][j] = 2;
			else if(dp[i][j] == (dp[i][j-1] + penalty))
				bt[i][j] = 3;
			else
				bt[i][j] = 4;
			
			if(dp[i][j] > result)
			{
				result = dp[i][j];
				query_end = i;
				database_end = j;
			}
		}
	}
	
	i = query_end;
	j = database_end;
	
	while(dp[i][j] > 0 && i > 0 && j > 0)
	{
		tempd = "";
		tempq = "";
		if(bt[i][j] == 1)
		{
			i--;
			j--;
			tempq.push_back(query[i]);
			tempd.push_back(database[j]);
			laq.insert(0,tempq);
			lad.insert(0,tempd);
		}
		else if(bt[i][j] == 2)
		{
			i--;
			tempq.push_back(query[i]);
			laq.insert(0,tempq);
			lad.insert(0,".");
		}
		else if(bt[i][j] == 3)
		{
			j--;
			tempd.push_back(database[j]);
			laq.insert(0,".");
			lad.insert(0,tempd);
		}
		else
			break;
	}	

	query_start = i+1;
	database_start = j+1;
	
	alignment res;
	
	res.query_id = query_id;
	res.database_id = database_id;
	res.query_sequence = laq;
	res.database_sequence = lad;
	res.query_start = query_start - 1;
	res.database_start = database_start - 1;
	res.query_end = query_end;
	res.database_end = database_end;
	res.score = result;
	
	for(int i = 0; i < query.length() + 1; ++i) {
        delete[] bt[i];   
    }
    
    delete[] bt;
	
	for(int i = 0; i < query.length() + 1; ++i) {
        delete[] dp[i];   
    }

    delete[] dp;
	
	return res;	
}

alignment global_alignment(int penalty,string query_id,string database_id, string query, string database,string alphabet,vector <vector <int> > &score_mat)
{
	int result = 0;
	
	string laq = "", lad = "",tempd = "",tempq = "";
	
	int **dp = 0;

    dp = new int *[query.length() + 1];
	
	for( int i = 0 ; i < query.length() + 1 ; i++ )
        dp[i] = new int [database.length() + 1];

	int **bt = 0;

    bt = new int *[query.length() + 1];
	
	for( int i = 0 ; i < query.length() + 1 ; i++ )
        bt[i] = new int [database.length() + 1];
	
	int i = 0, j = 0;

	for(i = 0;i <= query.length();i++)
		dp[i][0] = penalty*i;

	for(i = 0;i <= database.length();i++)
		dp[0][i] = penalty*i;

	for(i = 0;i <= query.length();i++)
		bt[i][0] = 2;

	for(i = 0;i <= database.length();i++)
		bt[0][i] = 3;

	for(int i = 1;i <= query.length(); i++)
	{
		for(int j = 1;j <= database.length();j++)
		{
			dp[i][j] = max(dp[i-1][j-1] + score(alphabet,score_mat,query[i-1],database[j-1]),max(dp[i-1][j] + penalty,dp[i][j-1] + penalty));
			
			if(dp[i][j] == (dp[i-1][j-1] + score(alphabet,score_mat,query[i-1],database[j-1])))
				bt[i][j] = 1;
			else if(dp[i][j] == (dp[i-1][j] + penalty))
				bt[i][j] = 2;
			else
				bt[i][j] = 3;
			
		}
	}
	
	result = dp[query.length()][database.length()];
	
	i = query.length();
	j = database.length();
	
	while(i > 0 && j > 0)
	{
		tempd = "";
		tempq = "";
		if(bt[i][j] == 1)
		{
			i--;
			j--;
			tempq.push_back(query[i]);
			tempd.push_back(database[j]);
			laq.insert(0,tempq);
			lad.insert(0,tempd);
		}
		else if(bt[i][j] == 2)
		{
			i--;
			tempq.push_back(query[i]);
			laq.insert(0,tempq);
			lad.insert(0,".");
		}
		else if(bt[i][j] == 3)
		{
			j--;
			tempd.push_back(database[j]);
			laq.insert(0,".");
			lad.insert(0,tempd);
		}
		else
			break;
	}	

	while(i > 0)
	{
		tempq = "";
		i--;
		tempq.push_back(query[i]);
		laq.insert(0,tempq);
		lad.insert(0,".");
	}
	
	while(j > 0)
	{
		tempd = "";
		j--;
		tempd.push_back(database[j]);
		laq.insert(0,".");
		lad.insert(0,tempd);
	}
	
	alignment res;
	
	res.query_id = query_id;
	res.database_id = database_id;
	res.query_sequence = laq;
	res.database_sequence = lad;
	res.query_start = 0;
	res.database_start = 0;
	res.query_end = query.length() - 1;
	res.database_end = database.length() - 1;
	res.score = result;
	
	for(int i = 0; i < query.length() + 1; ++i) {
        delete[] dp[i];   
    }
    
    delete[] dp;

	for(int i = 0; i < query.length() + 1; ++i) {
        delete[] bt[i];   
    }
    
    delete[] bt;

	return res;	
}

int main(int argc, char *argv[])
{
	
	int size = 1,i = 0,j = 0, input = atoi(argv[1]),count = atoi(argv[6]),penalty = atoi(argv[7]);
	vector< alignment > results;
	clock_t start,finish;
	alignment result;
	vector <vector <int> > score_mat;
	vector <int> temp;
	vector < vector <string> > database_lines;
	vector < vector <string> > query_lines;
	vector <string> temp_line;
	string filename;
	string line = "", queryline = "", databaseline = "",alphabet = "",str = "",id = "";
	ifstream infile (argv[2]);
	ofstream outfile("duplicate.txt");
	ifstream scoringmatrix (argv[5]);
	ifstream database_file (argv[3]);
	ifstream alphabet_file (argv[4]);
	
	if(alphabet_file.is_open())
	{
		getline (alphabet_file,alphabet);
	}
	alphabet_file.close();
	
	if(scoringmatrix.is_open())
	{
		while ( getline (scoringmatrix,line) )
		{
			j = 0;
			istringstream ss(line);
			string token;
			temp.clear();
				
			while (getline(ss, token,' '))
			{
				if(atoi(token.c_str()) != 0)
				{
					temp.push_back(atoi(token.c_str()));
				}
				
			}
			score_mat.push_back(temp);
		}
	}
	scoringmatrix.close();
	
	if (infile.is_open())
	{
		while ( getline (infile,line) )
		{
			istringstream ss(line);
			string token;
			temp_line.clear();
		
			if(line[0] == '>')
			{
				getline(ss, token,' ');
				outfile<<str;
				outfile.close();
				temp_line.push_back(id);
				temp_line.push_back(str);
				query_lines.push_back(temp_line);
				id = token.substr(5);
				str = "";
				filename = "query";
				filename += to_string(size);
				filename += ".txt";
				size++;
				outfile.open(filename.c_str());
				filename = "";
			}
			else
			{
				if(line[line.size() - 1] == '\r')
					line.erase(line.size()-1);
				str += line;
			}
		}
		infile.close();
		outfile<<str;
		outfile.close();
	
	}
	
	else cout << "Unable to open file"; 
	
	
	if (database_file.is_open())
	{
		while ( getline (database_file,line) )
		{
			istringstream ss(line);
			string token;
			temp_line.clear();
		
			if(line[0] == '>')
			{
				getline(ss, token,' ');
				temp_line.push_back(id);
				temp_line.push_back(databaseline);
				database_lines.push_back(temp_line);
				id = token.substr(5);
				databaseline = "";
			}
			else
			{
				if(line[line.size() - 1] == '\r')
					line.erase(line.size()-1);
				databaseline += line;
			}
		}
		
	}
	database_file.close();
	
	switch(input)
	{
		case 1:
			cout<<"Global Sequence Alignment:"<<endl;
			for (i = 1;i < query_lines.size();i++)
			{
				start = clock();
				for(j = 1;j < database_lines.size();j++)
				{
					temp.clear();
					result = global_alignment(penalty,query_lines[i][0],database_lines[j][0],query_lines[i][1],database_lines[j][1],alphabet,score_mat);
					results.push_back(result);
//					cout<<result.query_id<<"	"<<result.database_id<<"	"<<result.score<<endl;
				}
				finish = clock();
//				cout<<query_lines[i][0]<<"	"<<query_lines[i][1].size()<<"	"<<finish - start<<endl<<endl;
			}
			std::sort(results.begin(), results.end(), FirstColumnOnlyCmp());
			cout<<endl<<endl;
			for(i = 0;i < count;i++)
				cout<<"Score ="<<results[i].score<<endl<<results[i].query_id<<"	"<<results[i].query_start<<"	"<<results[i].query_sequence<<endl<<results[i].database_id<<"	"<<results[i].database_start<<"	"<<results[i].database_sequence<<endl<<endl;
		break;
		case 2:
			cout<<"Local Sequence Alignment:"<<endl;
			for (i = 1;i < query_lines.size();i++)
			{
				start = clock();
				for(j = 1;j < database_lines.size();j++)
				{
					temp.clear();
					result = local_alignment(penalty,query_lines[i][0],database_lines[j][0],query_lines[i][1],database_lines[j][1],alphabet,score_mat);
					results.push_back(result);
//					cout<<result.query_id<<"	"<<result.database_id<<"	"<<result.score<<endl;
				}
				finish = clock();
//				cout<<query_lines[i][0]<<"	"<<query_lines[i][1].size()<<"	"<<finish - start<<endl<<endl;
			}
			std::sort(results.begin(), results.end(), FirstColumnOnlyCmp());
			cout<<endl<<endl;
			for(i = 0;i < count;i++)
				cout<<"Score = "<<results[i].score<<endl<<results[i].query_id<<"	"<<results[i].query_start<<"	"<<results[i].query_sequence<<endl<<results[i].database_id<<"	"<<results[i].database_start<<"	"<<results[i].database_sequence<<endl<<endl;
		break;
		case 3: 
		cout<<"Dovetail Sequence Alignment:"<<endl;
			for (i = 1;i < query_lines.size();i++)
			{
				start = clock();
				for(j = 1;j < database_lines.size();j++)
				{
					temp.clear();
					result = dovetail_alignment(penalty,query_lines[i][0],database_lines[j][0],query_lines[i][1],database_lines[j][1],alphabet,score_mat);
					results.push_back(result);
//					cout<<result.query_id<<"	"<<result.database_id<<"	"<<result.score<<endl;
				}
				finish = clock();
//				cout<<query_lines[i][0]<<"	"<<query_lines[i][1].size()<<"	"<<finish - start<<endl<<endl;
			}
			std::sort(results.begin(), results.end(), FirstColumnOnlyCmp());
			cout<<endl<<endl;
			for(i = 0;i < count;i++)
				cout<<"Score = "<<results[i].score<<endl<<results[i].query_id<<"	"<<results[i].query_start<<"	"<<results[i].query_sequence<<endl<<results[i].database_id<<"	"<<results[i].database_start<<"	"<<results[i].database_sequence<<endl<<endl;
		break;
		default:
//			result = dovetail_alignment(-3,query_lines[1][0],database_lines[1][0],query_lines[1][1],database_lines[3][1],alphabet,score_mat);
//			cout<<result.query_id<<"	"<<result.database_id<<"	"<<result.score<<"	"<<result.query_start<<"	"<<result.query_end<<"	"<<result.database_start<<"	"<<result.database_end<<"	"<<endl<<result.query_sequence<<endl<<result.database_sequence<<endl<<endl;
			cout<<endl<<"Invalid Input";
	}
	cin>>size;
	return 0;
}
