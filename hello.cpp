#include <iostream>
#include <climits>
#include <chrono>
#include <fstream>
#include <cstring>
#include <cmath>
#include <cstdlib>
#include <cstdint>
#include <ctime>

using namespace std;

struct Stats 
{
    long long przypisania = 0;
    long long porownania  = 0; 
};

void fill_random(int* A, int n, unsigned int seed) {
    std::srand(seed);
    for (int i = 0; i < n; ++i)
        A[i] = (std::rand() % 2000001) - 1000000; // zakres [-10^6, 10^6]
}

//INSERTION SORT
void INSERTION_SORT(int A[], int n, Stats& stats)
{
    for (int i = 1; i < n; i++) {
        int x = A[i];                              
        int j = i - 1;
        while (j >= 0) {
            stats.porownania++;         
            if (A[j] > x) {
                A[j + 1] = A[j];        
                stats.przypisania++;
                j--;
            } else break;
        }
        A[j + 1] = x;
        stats.przypisania++;
    }
}

//INSERTION SORT modyfikacja
void INSERTION_SORT_DOUBLE(int A[], int n, Stats& stats)
{
    if (n < 2) return;

    //sortowanie pierwszej pary
    stats.porownania++;
    if (A[0] > A[1]) {
        int x = A[0];
        A[0] = A[1];  stats.przypisania++;
        A[1] = x;     stats.przypisania++;
    }

    //kolejne pary
    for (int i = 2; i < n; i += 2) {
        //gdy nie ma pelnej pary
        if (i + 1 >= n) {
            int x = A[i];
            int j = i - 1;
            while (j >= 0) {
                stats.porownania++;
                if (A[j] > x) {
                    A[j + 1] = A[j]; stats.przypisania++;
                    j--;
                } else break;
            }
            A[j + 1] = x; stats.przypisania++;
            break;
        }

        //nastepna para tak, by a <= b
        int a = A[i];     
        int b = A[i + 1]; 
        stats.porownania++;
        if (a > b) {
            int x = a; 
            a = b;     
            b = x;     
        }

        //wstaw większy (b) – przesuwamy o 2
        int j = i - 1;
        while (j >= 0) {
            stats.porownania++;
            if (A[j] > b) {
                A[j + 2] = A[j]; stats.przypisania++;
                j--;
            } else break;
        }
        A[j + 2] = b; stats.przypisania++;

        // wstaw mniejszy (a) – przesuwamy o 1
        while (j >= 0) {
            stats.porownania++;
            if (A[j] > a) {
                A[j + 1] = A[j]; stats.przypisania++;
                j--;
            } else break;
        }
        A[j + 1] = a; stats.przypisania++;
    }
}

//MERGE: 
void MERGE(int A[], int p, int s, int k, Stats& stats)
// jest dziwnie zapisane bo chcialem "swiadomie" zaznaczyc kiedy 
// uzywam danej zapisanej w stacku a kiedy w heapie (w sensie RAM, nie heapsort)
{
    int n1 = s - p + 1;
    int n2 = k - s;
 
    int* L = new int[n1 + 1]; // L jest wskaźnikiem do *L[0] na heapie
    int* R = new int[n2 + 1];

    // stack -> heap (kopiowanie z A do L/R)
    for (int i = 0; i < n1; ++i) { *(L + i) = A[p + i]; stats.przypisania++; }
    for (int j = 0; j < n2; ++j) { *(R + j) = A[s + 1 + j]; stats.przypisania++; }

    *(L + n1) = INT_MAX; stats.przypisania++;
    *(R + n2) = INT_MAX; stats.przypisania++;

    // scalanie: heap -> stack (z L/R do A)
    int i = 0, j = 0;
    for (int l = p; l <= k; ++l) {
        stats.porownania++; // porównanie *(L+i) <- *(R+j)
        if (*(L + i) <= *(R + j)) {
            A[l] = *(L + i); stats.przypisania++; // stack <- heap
            ++i;
        } else {
            A[l] = *(R + j); stats.przypisania++; // stack <- heap
            ++j;
        }
    }

    delete[] L;
    delete[] R;
}

//MERGE_SORT
void MERGE_SORT(int A[], int p, int k, Stats& stats)
{
    if (p >= k) return;
    int s = (p + k) / 2;
    // (nie liczymy porównań indeksów; brak przypisań danych)
    MERGE_SORT(A, p, s, stats);
    MERGE_SORT(A, s + 1, k, stats);
    MERGE(A, p, s, k, stats);
}

//MERGE_SORT modyfikacja
void MERGE_SORT3(int A[], int p, int k, Stats& stats)
{
    if (p >= k) return;
    int len = k - p + 1;
    int t1 = p + (len/3) - 1;
    int t2 = p + (2*len/3) - 1;

    // dopasowanie granic (gdy krótkie fragmenty)
    if (t1 < p) t1 = p;
    if (t2 < t1 + 1) t2 = t1 + 1;
    if (t2 > k - 1) t2 = k - 1;

    MERGE_SORT3(A, p,     t1, stats);
    MERGE_SORT3(A, t1+1,  t2, stats);
    MERGE_SORT3(A, t2+1,  k,  stats);

    MERGE(A, p,   t1, t2, stats); // połączenie 1
    MERGE(A, p,   t2, k,  stats); // połączenie 2
}

//HEAP SORT
void swap_count(int& x, int& y, Stats& stats) {
    int tmp = x;
    x = y;       stats.przypisania++;
    y = tmp;     stats.przypisania++;
}

void heapify_binary(int A[], int n, int i, Stats& stats)
{
    int l = 2*i + 1;
    int r = 2*i + 2;
    int largest = i;

    if (l < n) { stats.porownania++; if (A[l] > A[largest]) largest = l; }
    if (r < n) { stats.porownania++; if (A[r] > A[largest]) largest = r; }

    if (largest != i) {
        swap_count(A[i], A[largest], stats);
        heapify_binary(A, n, largest, stats);
    }
}

void HEAP_SORT(int A[], int n, Stats& stats)
{
    // build-heap (ostatni rodzic = n/2 - 1)
    for (int i = n/2 - 1; i >= 0; --i)
        heapify_binary(A, n, i, stats);

    // sortowanie
    for (int i = n - 1; i > 0; --i) {
        swap_count(A[0], A[i], stats);
        heapify_binary(A, i, 0, stats);
    }
}

//HEAP SORT modyfikacja

void heapify_ternary(int A[], int n, int i, Stats& stats)
{
    int c1 = 3*i + 1;
    int c2 = 3*i + 2;
    int c3 = 3*i + 3;
    int largest = i;

    if (c1 < n) { stats.porownania++; if (A[c1] > A[largest]) largest = c1; }
    if (c2 < n) { stats.porownania++; if (A[c2] > A[largest]) largest = c2; }
    if (c3 < n) { stats.porownania++; if (A[c3] > A[largest]) largest = c3; }

    if (largest != i) {
        swap_count(A[i], A[largest], stats);
        heapify_ternary(A, n, largest, stats);
    }
}

void HEAP_SORT_TERNARY(int A[], int n, Stats& stats)
{
    for (int i = (n - 2) / 3; i >= 0; --i)
        heapify_ternary(A, n, i, stats);

    for (int i = n - 1; i > 0; --i) {
        swap_count(A[0], A[i], stats);
        heapify_ternary(A, i, 0, stats);
    }
}

// funkcja testujaca
/*void print_array(const char* title, const int A[], int n) {
    cout << title;
    for (int i = 0; i < n; i++) cout << A[i] << (i+1<n?' ':'\n');
}*/

int main() {
    const int Ns[] = {100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000};
    const int R = 10;
    std::srand(std::time(NULL)); // dla kazdego wywolania programu inny seed
    unsigned int base_seed = std::rand(); // na podstawie tego powyzej, ktory jest losowany na podstawie czasu
                                            // mamy base_seed ktorego uzywamy juz do kazdej z dziesieciu probs

     // CSV średnie 
    ofstream out("wyniki.csv", ios::trunc);
    out << "n,"
        << "INS_czas_ns,INS_por,INS_przyp,"
        << "INSD_czas_ns,INSD_por,INSD_przyp,"
        << "M2_czas_ns,M2_por,M2_przyp,"
        << "M3_czas_ns,M3_por,M3_przyp,"
        << "H2_czas_ns,H2_por,H2_przyp,"
        << "H3_czas_ns,H3_por,H3_przyp\n";




        for (int n : Ns) {
        long double sum_tINS=0, sum_cINS=0, sum_aINS=0; // dla kazdego n z Ns zerujemy najpierw: dla kazdego algorytmu:
        long double sum_tINSD=0, sum_cINSD=0, sum_aINSD=0;                      // czasy, porownania, przypisania
        long double sum_tM2=0, sum_cM2=0, sum_aM2=0;
        long double sum_tM3=0, sum_cM3=0, sum_aM3=0;
        long double sum_tH2=0, sum_cH2=0, sum_aH2=0;
        long double sum_tH3=0, sum_cH3=0, sum_aH3=0;

        for (int r = 0; r < R; ++r) {
            int* baza = new int[n];
            fill_random(baza, n, base_seed + n*37 + r);

           
            // 1) INSERTION_SORT
            
            int* A1 = new int[n];
            memcpy(A1, baza, n*sizeof(int)); //kopiujemy baze losowych do A1 - kazdy algorytm ma swoja A_n ktora sortuje - w ten sposob nie dotyka bazy 
            Stats sINS;
            auto t1 = chrono::high_resolution_clock::now();
            INSERTION_SORT(A1, n, sINS);
            auto t2 = chrono::high_resolution_clock::now();
            long long tINS = chrono::duration_cast<chrono::nanoseconds>(t2-t1).count();
            delete[] A1;

            
            // 2) INSERTION_SORT_DOUBLE
            
            int* A2 = new int[n];
            memcpy(A2, baza, n*sizeof(int)); 
            Stats sINSD;
            t1 = chrono::high_resolution_clock::now();
            INSERTION_SORT_DOUBLE(A2, n, sINSD);
            t2 = chrono::high_resolution_clock::now();
            long long tINSD = chrono::duration_cast<chrono::nanoseconds>(t2-t1).count();
            delete[] A2;

            
            // 3) MERGE_SORT
            
            int* A3 = new int[n];
            memcpy(A3, baza, n*sizeof(int));
            Stats sM2;
            t1 = chrono::high_resolution_clock::now();
            MERGE_SORT(A3, 0, n-1, sM2);
            t2 = chrono::high_resolution_clock::now();
            long long tM2 = chrono::duration_cast<chrono::nanoseconds>(t2-t1).count();
            delete[] A3;

            
            // 4) MERGE_SORT3
            
            int* A4 = new int[n];
            memcpy(A4, baza, n*sizeof(int));
            Stats sM3;
            t1 = chrono::high_resolution_clock::now();
            MERGE_SORT3(A4, 0, n-1, sM3);
            t2 = chrono::high_resolution_clock::now();
            long long tM3 = chrono::duration_cast<chrono::nanoseconds>(t2-t1).count();
            delete[] A4;

            // 5) HEAP_SORT (binarny)
            
            int* A5 = new int[n];
            memcpy(A5, baza, n*sizeof(int));
            Stats sH2;
            t1 = chrono::high_resolution_clock::now();
            HEAP_SORT(A5, n, sH2);
            t2 = chrono::high_resolution_clock::now();
            long long tH2 = chrono::duration_cast<chrono::nanoseconds>(t2-t1).count();
            delete[] A5;

            
            // 6) HEAP_SORT_TERNARY
            
            int* A6 = new int[n];
            memcpy(A6, baza, n*sizeof(int));
            Stats sH3;
            t1 = chrono::high_resolution_clock::now();
            HEAP_SORT_TERNARY(A6, n, sH3);
            t2 = chrono::high_resolution_clock::now();
            long long tH3 = chrono::duration_cast<chrono::nanoseconds>(t2-t1).count();
            delete[] A6;

            // kazde z 10 powtorzen dla ustalonego n dodajemy: czasy, porownania przypisania dla kazdego z algorytmow
            sum_tINS+=tINS;   sum_cINS+=sINS.porownania;   sum_aINS+=sINS.przypisania;
            sum_tINSD+=tINSD; sum_cINSD+=sINSD.porownania; sum_aINSD+=sINSD.przypisania;
            sum_tM2+=tM2;     sum_cM2+=sM2.porownania;     sum_aM2+=sM2.przypisania;
            sum_tM3+=tM3;     sum_cM3+=sM3.porownania;     sum_aM3+=sM3.przypisania;
            sum_tH2+=tH2;     sum_cH2+=sH2.porownania;     sum_aH2+=sH2.przypisania;
            sum_tH3+=tH3;     sum_cH3+=sH3.porownania;     sum_aH3+=sH3.przypisania;

            delete[] baza;
        }

        auto avg = [R](long double x)
        { return (long long) llround(x / R); }; // lambda przyjmuje long double x=sum_tALGORYTM, i uzywa [R] z pętli
                                                    // ma dostep do R bo ta zmienna jest w main() i lambda tez jest w main()
        out << n << ","
            << avg(sum_tINS)  << "," << avg(sum_cINS)  << "," << avg(sum_aINS)  << ","
            << avg(sum_tINSD) << "," << avg(sum_cINSD) << "," << avg(sum_aINSD) << ","
            << avg(sum_tM2)   << "," << avg(sum_cM2)   << "," << avg(sum_aM2)   << ","
            << avg(sum_tM3)   << "," << avg(sum_cM3)   << "," << avg(sum_aM3)   << ","
            << avg(sum_tH2)   << "," << avg(sum_cH2)   << "," << avg(sum_aH2)   << ","
            << avg(sum_tH3)   << "," << avg(sum_cH3)   << "," << avg(sum_aH3)
            << "\n"; // kazde wejscie w jednym wierszu to po prostu suma wszystkich (czasów/porownan/przypisan) podzielona przez 10 z uzyciem llround

        cout << "n=" << n << " zakończono.\n";
    }

    cout << "Zapisano wyniki do wyniki.csv\n";


    /*
    // MERGE_SORT
    int B[] = {5,3,4,1,2};
    int b = sizeof(B)/sizeof(B[0]);
    Stats merge_stats;
    MERGE_SORT(B, 0, b-1, merge_stats);
    print_array("MERGE_SORT: ", B, b);
    cout << "MERGE_SORT stats -> porownania=" << merge_stats.porownania
         << ", przypisania=" << merge_stats.przypisania << "\n\n";

    // MERGE_SORT3
    int C[] = {9,1,8,2,7,3,6,4,5};
    int c = sizeof(C)/sizeof(C[0]);
    Stats merge3_stats;
    MERGE_SORT3(C, 0, c-1, merge3_stats);
    print_array("MERGE_SORT3: ", C, c);
    cout << "MERGE_SORT3 stats -> porownania=" << merge3_stats.porownania
         << ", przypisania=" << merge3_stats.przypisania << "\n\n";

    // INSERTION_SORT
    int A1[] = {2, 35, 4, 6};
    int n1 = sizeof(A1)/sizeof(A1[0]);
    Stats ins_stats;
    INSERTION_SORT(A1, n1, ins_stats);
    print_array("INSERTION_SORT: ", A1, n1);
    cout << "INSERTION_SORT stats -> porownania=" << ins_stats.porownania
         << ", przypisania=" << ins_stats.przypisania << "\n\n";

    // INSERTION_SORT_DOUBLE (modyfikacja)
    int A2[] = {7,3,5,1,9,0,8};
    int n2 = sizeof(A2)/sizeof(A2[0]);
    Stats ins2_stats;
    INSERTION_SORT_DOUBLE(A2, n2, ins2_stats);
    print_array("INSERTION_SORT_DOUBLE: ", A2, n2);
    cout << "INSERTION_SORT_DOUBLE stats -> porownania=" << ins2_stats.porownania
         << ", przypisania=" << ins2_stats.przypisania << "\n\n";

    // HEAP_SORT (binarny)
    int H1[] = {10,5,3,2,4,7,8,9,1,6};
    int h1 = sizeof(H1)/sizeof(H1[0]);
    Stats heap_stats;
    HEAP_SORT(H1, h1, heap_stats);
    print_array("HEAP_SORT (binarny): ", H1, h1);
    cout << "HEAP_SORT stats -> porownania=" << heap_stats.porownania
         << ", przypisania=" << heap_stats.przypisania << "\n\n";

    // HEAP_SORT_TERNARY modyfikacja
    int H2[] = {10,5,3,2,4,7,8,9,1,6,11,0,12};
    int h2 = sizeof(H2)/sizeof(H2[0]);
    Stats heap3_stats;
    HEAP_SORT_TERNARY(H2, h2, heap3_stats);
    print_array("HEAP_SORT (ternarny): ", H2, h2);
    cout << "HEAP_SORT_TERNARY stats -> porownania=" << heap3_stats.porownania
         << ", przypisania=" << heap3_stats.przypisania << "\n";
    */

    return 0;
}
