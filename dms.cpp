#include<bits/stdc++.h>
using namespace std;

//四元组 存储了一个r_mer以及它所在的序列号，在序列中的位置，以及类型（0表示r_mer 1表示q_mer的d_neighbor） 
class Rmer
{
	public:
		string str;
		int seq_number, pos, type;
		
		Rmer(string s, int seqnum, int position, int typii): seq_number(seqnum), str(s), pos(position), type(typii) {}
		
};

const string ATGC="ATGC";


/*  重要的参数  */ 

//模体长度 编辑距离 
int p=10, d=2;
//生成数据的序列数目和序列长度  由于bitset的原因，序列数目需要是const的 
const int num=5, length=600;

/*  参数结束  */ 


//数据 
vector<string> data;
// 
vector<Rmer> p_nbr, r_mer;
//统计矩阵 
vector<vector<bitset<num+1> > > A;

//随机生成数据，存入data中，数据的参数参见上面注释 
//可以控制保存和加载数据
//数据文件data.txt与代码文件在同一目录 
void data_generation(bool save=true, bool load=false)
{
	string str;
	if(load)
	{
		fstream fin("data.txt", ios::in);
		for(int i=0; i<num; i++)
		{
			fin >> str;
			data.push_back(str);
		}
		fin.close();
	}
	else
	{
		for(int i=0; i<num; i++)
		{
			str = "";
			for(int j=0; j<length; j++)
			    str += ATGC[rand()%4];
			data.push_back(str);
		}
	}
	if(save)
	{
		fstream fout("data.txt", ios::out);
		for(int i=0; i<num; i++)
			fout << data[i] << endl;
		fout.close();
	}
}

// 动态规划过程 
// 原本用L[i][j]表示串V的前i个字符和串W的前j个字符之间最小的编辑距离
// 现在这里又滚动数组的意思 现在一维的L表示原来的L[i-1]，而S表示L[i]， 
vector<int> row(vector<int>&L, char& a, string& str)
{
	vector<int> S;
	S.push_back(L[0]+1);
	for(int i=1; i<L.size(); i++)
		//这几个代表哪几种情况来着？  V加_  W加_  对上 
		S.push_back(min(S[i-1]+1, min(L[i]+1, L[i-1]+(a==str[i-1]? 0:1))));
	return S;
}

void gen(string& p_mer, int i, vector<int> L, string& V, int& seqnum, int& position)
{
	for(auto a:ATGC)
	{
	    auto Li = row(L, a, p_mer);
		V[i-1] = a;
		if(Li[p] <= d)
		    p_nbr.push_back(Rmer(V.substr(0, i), seqnum, position, 1));
//		else
//		{
			// 外面那圈else似乎没有意义 
			
			if(*min_element(Li.begin(), Li.end()) <= d)
			    gen(p_mer, i+1, Li, V, seqnum, position);
//		},
	}	
}

void word_neighbor_generation(string& p_mer, int& seqnum, int& position)
{
	vector<int> L;
	string V;
	V.resize(p+d+2);
	for(int i=0; i<=p; i++)
	    L.push_back(i);
	gen(p_mer, 1, L, V, seqnum, position);
}

void p_mer_neighbor_generation()
{
	for(int i=0; i<data.size(); i++)
	{
		for(int j=0; j+p<=data[i].length(); j++)
		{
			string s = data[i].substr(j, p);
			//cout << s<< endl;
			// 貌似传引用的时候不能传一个临时产生的串 
			word_neighbor_generation(s, i, j);
		}
		    
	}
}

void r_mer_generation()
{
	for(int i=0; i<data.size(); i++)
	{
		for(int j=0; j+p-d<=data[i].length(); j++)
		    for(int k=p-d; k<=p+d&&j+k<=data[i].length(); k++)
		        r_mer.push_back(Rmer(data[i].substr(j, k), i, j, 0));
	}
}

void do_dms()
{
	// 按照前两项对r_mer进行排序 
	sort(r_mer.begin(), r_mer.end(), [](const Rmer&a, const Rmer&b) {if(a.str<b.str) return true; else if(a.str==b.str&&a.seq_number<b.seq_number) return true; else return false;});
	auto it=r_mer.begin();
	// 删除前两项重复的r_mer，因为对于一个r_mer只需要记录它在某条序列里出现过就可以了 
	while(it != r_mer.end())
	{
		auto left=it+1, right=it+1;
		while(right!=r_mer.end()&&right->str==it->str&&right->seq_number==it->seq_number)
		    right++;
		it = r_mer.erase(left, right);
	}
	// 合并r_mer和所以p_mer的d_neighbor 
	r_mer.insert(r_mer.end(), p_nbr.begin(), p_nbr.end());
	// 合并后按照第一个和最后一个关键字进行排序 
	sort(r_mer.begin(), r_mer.end(), [](const Rmer&a, const Rmer&b) {if(a.str<b.str) return true; else if(a.str==b.str&&a.type<b.type) return true; else return false;});
	A.resize(num+1);
	for(int i=0; i<A.size(); i++)
	    A[i].resize(length);
	// 计算过程
	// A[i][j][k]=1表示序列编号i的序列中起始位置为j的p_mer在序列编号为k的序列中存在一个d_neighbor 
	for(int i=0; i<r_mer.size(); )
	{
		int j=i, sep=-1;
		while(j<r_mer.size()&&r_mer[j].str==r_mer[i].str)
		{
			if(r_mer[j].type==0) sep = j;
			j++;
		}
		if(sep!=-1)
		{
			for(int k=sep+1; k<j; k++)
		        for(int l=i; l<=sep; l++)
		        {
//                  这两句注释是我验证结果用的，输出某个模体在所有r_mer中的出现情况		        	
//		    	    if(r_mer[k].seq_number==0&&r_mer[k].pos==11)
//		    	        cout << r_mer[l].seq_number << " " << r_mer[l].str << ' ' << r_mer[l].type << endl;
				    A[r_mer[k].seq_number][r_mer[k].pos][r_mer[l].seq_number] = 1;
			    }
		}      
		i = j;
	}
	// 找到的模体存储在q_mer.txt中，也跟代码在同一目录下
	// 输出的一行表示一个q_mer，一行两个数字，分别表示序列编号和起始位置 
	fstream fout("q_mer.txt", ios::out);
	// 如果对于一个p_mer在每条序列中都能找到它的d_neighbor，则为一个模体 
	for(int i=0; i<A.size(); i++)
	    for(int j=0; j<A[i].size(); j++)
	        if(A[i][j].count() == num)
	            fout << i << '\t' << j << endl;
	fout.close();
}

int main()
{
	data_generation(); //false, true    要对同一组生成的数据进行重复测试时把这两个参数传进去，实际上也可以直接注释掉时间做种 
	cout << "data got\n";
	r_mer_generation();
	cout << "r_mer got\n";
	p_mer_neighbor_generation();
	cout << "p_mer_d_neighbor got\n";
	do_dms();
	cout << "done\n";
	return 0;
}
