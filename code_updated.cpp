#include <bits/stdc++.h>
using namespace std;

#define int long long
#define uid(a, b) uniform_int_distribution<int>(a, b)(mt)

int cmps=0;
int failed=0;

//partition when we have just one pivot. used for quickselect.
int one_pivot_partition(vector<int> &a, int l, int r, pair<int,int> mid)
{
    vector<int> less,more;
    int pvt_idx;
    for(int i=l;i<=r;i++)
    {
        if(i==mid.second)
            continue;
        cmps++;
        if(a[i]<mid.first)
            less.push_back(a[i]);
        else
            more.push_back(a[i]);
    }
    int ptr=l;
    for(int i:less)
        a[ptr]=i, ptr++;
    pvt_idx=ptr;
    a[ptr]=mid.first;
    ptr++;
    for(int i:more)
        a[ptr]=i, ptr++;
    return pvt_idx;
}

//partition function. takes an array a, range [l, r], two pivots [a[i1], i1], [a[i2], i2] and partitions the array based on these pivots, and returns the new position (that is, ranks) of the pivots in subarray a[l:r]
pair<int,int> partition(vector<int> &a, int l, int r, pair<int,int> mid1, pair<int,int> mid2)
{
    if(mid1==mid2)
    {
        int pvt_idx=one_pivot_partition(a, l, r, mid1);
        return make_pair(pvt_idx, pvt_idx);
    }
    pair<int, int> new_pivots;
    vector<int> less,more, bw;//less->elements in a less than smaller pivot. similarly bw->between and more->more than larger pivot
    int smaller_idx=min(mid1.second, mid2.second), larger_idx=max(mid1.second, mid2.second);
    for(int i=l;i<=r;i++)
    {
        if(i==smaller_idx || i==larger_idx)//these are not counted in comparisons, as they are not comparing array elements. Also, these comparisons can be easily avoided
            continue;
        if(a[i]<mid1.first)
            cmps+=1, less.push_back(a[i]);
        else if(a[i]>=mid2.first)
            cmps+=2, more.push_back(a[i]);
        else
            cmps+=2, bw.push_back(a[i]);
    }
    int ptr=l;
    for(int i:less)
        a[ptr]=i, ptr++;
    new_pivots.first=ptr;
    a[ptr]=mid1.first;
    ptr++;
    for(int i:bw)
        a[ptr]=i, ptr++;
    new_pivots.second=ptr;
    a[ptr]=mid2.first;
    ptr++;
    for(int i:more)
        a[ptr]=i, ptr++;
    return new_pivots;
}


//generates random pivot for quicksort function. uses randomization to ensure that pivot is selected near median for faster convergence
pair<int,int> gen_pivot(vector<int> &a, int l, int r)
{
    if(l==r)
        return make_pair(a[l],l);
    random_device rd;
    mt19937 mt(rd());
    int n = r - l + 1;
    int k= min(n/2, (int)(10*log2(n)));
    vector<pair<int, int>> s;
    for(int i=0;i<k;i++)
    {
        int idx = uid(0, n - 1);
        s.push_back(make_pair(a[idx + l], idx+l));
    }
    sort(s.begin(),s.end());
    cmps+=k*ceil(log2(k));//adding comparisons of merge sort, as it takes exactly k*ceil(log2(k)) comparisons
    return s[(k-1)/2];
}

int quickSelect(int l, int r, vector<int>&a, int k)//finds kth element in a[l:r]. l<=k<=r
{
    if(l < r)
    {
        pair<int,int> mid_pvt=gen_pivot(a, l, r);//generates partition. also ensures that pivot lies near middle as an optimization.
        int pvt_idx = one_pivot_partition(a, l, r, mid_pvt);//uses single pivot partition
        // int pvt_idx = partition(a, l, r, mid_pvt, mid_pvt).first;//uses single pivot partition
        cmps++;
        if (pvt_idx == k)
            return a[pvt_idx];
        cmps++;
        if(pvt_idx > k)
            return quickSelect(l, pvt_idx-1, a, k);
        else
            return quickSelect(pvt_idx+1, r, a, k);
    }
    else if(l==r)
        return a[k];
    return 0;
}

//returns the subarray C (from doc) and index of median in that subarray.
pair<int, vector<int>> reduce(vector<int> &a, float alpha, float beta)//n^alpha=subarray size that we desire. n^beta =number of elements seleccted to select pivots using randomization
{
    int n=a.size();
    random_device rd;
    mt19937 mt(rd());
    int k=pow(n, beta);
    vector<int> s;
    for(int i=0;i<k;i++)
    {
        int idx = uid(0, n - 1);
        s.push_back(((a[idx])<<32)+idx);//this is done to encode index so that quickSelect and then partition can be used
    }
    int idx1=max(0ll, (int)(  (k-1)/2-k/(2*pow(n, 1-alpha))   )    );//rank of first pivot in set s
    int idx2=min(k-1ll, (int)(   (k-1)/2+k/(2*pow(n, 1-alpha))  )   );//rank of second pivot in set s
    int pvt1=quickSelect(0, k-1, s, idx1);//first pivot
    int pvt2=quickSelect(0, k-1, s, idx2);//second pivot

    pair<int,int> pivots=partition(a, 0, n-1, make_pair(pvt1>>32, pvt1%(1ll<<32)), make_pair(pvt2>>32, pvt2%(1ll<<32)));//partitioning the array based on pivots
    vector<int> reduced;//subarray of o(n) size  which can be used to get median.
    for(int i=pivots.first;i<=pivots.second;i++)
    {
        reduced.push_back(a[i]);
    }
    int idx=(n-1)/2-pivots.first;
    return make_pair(idx, reduced);
}

int median(vector<int> &a, float alpha, float beta)//the main function for calculating the median
{
    int n=a.size();
    auto ret=reduce(a, alpha, beta);//reducing the array to o(n) size
    vector<int> mid_array=ret.second;
    int idx=ret.first;
    int sz=mid_array.size();
    if(idx<0 || idx>=mid_array.size())//if median does not lie in the subarray, use quickselect (or sorting) to find median
    {
        failed++;
        return quickSelect(0, n-1, a, (n-1)/2);
    }
    return quickSelect(0, sz-1, mid_array, idx);
}



void func(float alpha, float beta)
{
    random_device rd;
    mt19937 gen(rd());
    int n=30000000;
    vector<int> b(n);
    for(int i=0;i<n;i++)
        b[i]=i;

    shuffle(b.begin(), b.end(), gen);//We have used shuffle function to generate random permutation
    float fin_ans=0;
    vector<int> a=b;
    sort(a.begin(),a.end());//just for finding the true median to check if our median is correct
    int med=a[(n-1)/2];//actual median
    for(int i=0;i<100;i++)
    {
        cmps=0;
        failed=0;
        vector<int> a=b;//a again set to orignal array
        int my_med=median(a, alpha, beta);
        if(my_med!=med)//checks if median we calculated is correct or nor
        {
            cout<<"Incorrect median calculated for iteration i="<<i<<", actual median: "<<med<<", your algo's median: "<<my_med<<endl;
            exit(-1);
        }
        fin_ans+=(float)cmps/n;
        cout<<(float)cmps/n<<endl;
    }
    cout<<"a="<<alpha<<endl<<"b="<<beta<<endl<<"num_cmps/n="<<fin_ans/100<<endl;
    cout<<"failed_cnt: "<<failed<<endl<<endl;
}
//Running the code will print average number of comparisons divided by n for n=30000000. a and b can be changed from function call from main function. n can be changed from above. Also, number of iterations can be changed.
//failed count will give number of times (from the 50 iterations) when median does not lie in selected subarray (and thus quickselect is used in that case)

//The code may take upto a minute to execute, as array size and number of iterations are huge.



//NOTE: In code, to calculate ith ranked element for o(n) sized arrays (the set B and subarray C in doc) quickSelect is used instead of sorting as an implementation optimization.
signed main()
{
    func(0.7,0.7);
    return 0;
}