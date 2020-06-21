#include<bits/stdc++.h>
using namespace std;

//��Ԫ�� �洢��һ��r_mer�Լ������ڵ����кţ��������е�λ�ã��Լ����ͣ�0��ʾr_mer 1��ʾq_mer��d_neighbor�� 
class Rmer
{
	public:
		string str;
		int seq_number, pos, type;
		
		Rmer(string s, int seqnum, int position, int typii): seq_number(seqnum), str(s), pos(position), type(typii) {}
		
};

const string ATGC="ATGC";


/*  ��Ҫ�Ĳ���  */ 

//ģ�峤�� �༭���� 
int p=10, d=2;
//�������ݵ�������Ŀ�����г���  ����bitset��ԭ��������Ŀ��Ҫ��const�� 
const int num=5, length=600;

/*  ��������  */ 


//���� 
vector<string> data;
// 
vector<Rmer> p_nbr, r_mer;
//ͳ�ƾ��� 
vector<vector<bitset<num+1> > > A;

//����������ݣ�����data�У����ݵĲ����μ�����ע�� 
//���Կ��Ʊ���ͼ�������
//�����ļ�data.txt������ļ���ͬһĿ¼ 
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

// ��̬�滮���� 
// ԭ����L[i][j]��ʾ��V��ǰi���ַ��ʹ�W��ǰj���ַ�֮����С�ı༭����
// ���������ֹ����������˼ ����һά��L��ʾԭ����L[i-1]����S��ʾL[i]�� 
vector<int> row(vector<int>&L, char& a, string& str)
{
	vector<int> S;
	S.push_back(L[0]+1);
	for(int i=1; i<L.size(); i++)
		//�⼸�������ļ���������ţ�  V��_  W��_  ���� 
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
			// ������Ȧelse�ƺ�û������ 
			
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
			// ò�ƴ����õ�ʱ���ܴ�һ����ʱ�����Ĵ� 
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
	// ����ǰ�����r_mer�������� 
	sort(r_mer.begin(), r_mer.end(), [](const Rmer&a, const Rmer&b) {if(a.str<b.str) return true; else if(a.str==b.str&&a.seq_number<b.seq_number) return true; else return false;});
	auto it=r_mer.begin();
	// ɾ��ǰ�����ظ���r_mer����Ϊ����һ��r_merֻ��Ҫ��¼����ĳ����������ֹ��Ϳ����� 
	while(it != r_mer.end())
	{
		auto left=it+1, right=it+1;
		while(right!=r_mer.end()&&right->str==it->str&&right->seq_number==it->seq_number)
		    right++;
		it = r_mer.erase(left, right);
	}
	// �ϲ�r_mer������p_mer��d_neighbor 
	r_mer.insert(r_mer.end(), p_nbr.begin(), p_nbr.end());
	// �ϲ����յ�һ�������һ���ؼ��ֽ������� 
	sort(r_mer.begin(), r_mer.end(), [](const Rmer&a, const Rmer&b) {if(a.str<b.str) return true; else if(a.str==b.str&&a.type<b.type) return true; else return false;});
	A.resize(num+1);
	for(int i=0; i<A.size(); i++)
	    A[i].resize(length);
	// �������
	// A[i][j][k]=1��ʾ���б��i����������ʼλ��Ϊj��p_mer�����б��Ϊk�������д���һ��d_neighbor 
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
//                  ������ע��������֤����õģ����ĳ��ģ��������r_mer�еĳ������		        	
//		    	    if(r_mer[k].seq_number==0&&r_mer[k].pos==11)
//		    	        cout << r_mer[l].seq_number << " " << r_mer[l].str << ' ' << r_mer[l].type << endl;
				    A[r_mer[k].seq_number][r_mer[k].pos][r_mer[l].seq_number] = 1;
			    }
		}      
		i = j;
	}
	// �ҵ���ģ��洢��q_mer.txt�У�Ҳ��������ͬһĿ¼��
	// �����һ�б�ʾһ��q_mer��һ���������֣��ֱ��ʾ���б�ź���ʼλ�� 
	fstream fout("q_mer.txt", ios::out);
	// �������һ��p_mer��ÿ�������ж����ҵ�����d_neighbor����Ϊһ��ģ�� 
	for(int i=0; i<A.size(); i++)
	    for(int j=0; j<A[i].size(); j++)
	        if(A[i][j].count() == num)
	            fout << i << '\t' << j << endl;
	fout.close();
}

int main()
{
	data_generation(); //false, true    Ҫ��ͬһ�����ɵ����ݽ����ظ�����ʱ����������������ȥ��ʵ����Ҳ����ֱ��ע�͵�ʱ������ 
	cout << "data got\n";
	r_mer_generation();
	cout << "r_mer got\n";
	p_mer_neighbor_generation();
	cout << "p_mer_d_neighbor got\n";
	do_dms();
	cout << "done\n";
	return 0;
}
