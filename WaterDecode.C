#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>

#define MAX_N 400
#define MAX_Nr 430
#define MAX_D 330

long double *F_store;
long double *B_store;

// Define a struct to hold necessary data
typedef struct {
    double Pi, Pd, Pt;
    double Ps;
    int N, Nr, Nd;
    int *info;
    int *r, *w, *d;
    double *priors;
    int *data;
} ChannelData;

long double forward(ChannelData *data, int i, int j);
long double backward(ChannelData *data, int i, int j);

void *chk_alloc
( unsigned n,		/* Number of elements */
  unsigned size		/* Size of each element */
)
{ 
  void *p;

  p = calloc(n,size);

  if (p==0)
  { fprintf(stderr,"Ran out of memory (while trying to allocate %d bytes)\n",
      n*size);
    exit(1);
  }

  return p;
}

bool over_limit(int i, int j){
    int a  = abs(i - j);
    if (a < 10){
        return false;
    }
    float max_drift = pow(3 * i, 0.5) / 2.1;
    if (a > max_drift){
        return true;
    }
    else{
        return false;
    }
}

void initialise_channel_data(ChannelData *data) {
    int s = (data->N + 1) * (data->Nr + 1);
    for (int i = 0; i < s; i++) {
        F_store[i] = -1.0;
        B_store[i] = -1.0;
    }
}

int *int_read_from_file(const char *filename, int size) {
    // fprintf(stderr, "Reading from file %s size %d \n", filename, size);
    FILE *file = fopen(filename, "r");
    if (!file) {
        printf("Error opening file: %s\n", filename);
        exit(1);
    }
    int *array = (int *)malloc(size * sizeof(int));
    int c;
    for (int i = 0; i < size; i++) {
        c = fscanf(file, "%d", &array[i]);
    if (c==EOF) 
      { if (i>0)
        { fprintf(stderr,
          "Warning: Short block (%d long) in %s \n",i, &filename);
        }
      }
    }

    fclose(file);
    return array;
}

double *d_read_from_file(const char *filename, int size) {
    // fprintf(stderr, "Reading from file %s with size %d \n", filename, size);
    FILE *file = fopen(filename, "r");
    if (!file) {
        printf("Error opening file: %s\n", filename);
        exit(1);
    }
    int c;
    double *array = (double *)malloc(size * sizeof(double));

    for (int i = 0; i < size; i++) {
        fscanf(file, "%lf", &array[i]);
    }
    return array;

    // for(int i = 0; i<30; i++) {
    //     fprintf(stderr, "%lf", array[i]);
    // }
    
    // for (int i = 0; i < size; i++) {
    //     c = fscanf(file, "%f", &my_arr[i]);
    //     fprintf(stderr, "%d", c);
    //     fprintf(stderr, "Found %g \n", my_arr[i]);
    //         if (c==EOF) 
    //   { if (i>0)
    //     { fprintf(stderr,
    //       "Warning: Short block (%d long) in %c\n",i, &filename);
    //     }
    //   }
    // }
    // fclose(file);
}


double gamm(ChannelData *data, int i1, int j1, int i2, int j2) {
    double g = 0.0;
    if (i2 == i1 && j2 == j1 + 1) {
        g = data->Pi / 2;
    } else if (i2 == i1 + 1 && j2 == j1) {
        g = data->Pd;
    } else if (i2 == i1 + 1 && j2 == j1 + 1) {
        bool found = false;
        int ind;
            for (int k=0; k<data->Nd; k++) 
            { 
            if (data->info[k] == i1) 
            { 
                found = true;
                ind = k;
                break; 
            } 
            } 
        if (!found) {
            g = (data->r[j1] == data->w[i1]) ? data->Pt * (1 - data->Ps) : data->Pt * data->Ps;
        } else {
            double p = (1.0 - data->priors[ind]); //probability that bit is 0
            g = ((data->r[j1] + data->w[i1]) % 2 == 1) ? 
                (p * data->Pt * data->Ps + (1 - p) * data->Pt * (1 - data->Ps)) :
                (p * data->Pt * (1 - data->Ps) + (1 - p) * data->Pt * data->Ps);
        }
    }
    if (g > 1){
        fprintf(stderr, "Got large gamma value");
    }
    return g;
}

long double compute_forward(ChannelData *data, int i, int j) {
    // fprintf(stderr, "Computing forwards for %d %d \n", i, j);
    long double f;
    if (i > data->N || j > data->Nr || i < 0 || j < 0) {
        fprintf(stderr, "Outside\n" );
        f = 0.0;
    } else if (i == 0 && j == 0) {
        // fprintf(stderr, "Boundary\n" );
        f  =  1.0;
        // fprintf(stderr, "Got forward value of %lf for %d %d \n", f, i, j);
    } else if (i == 0) {
        f = forward(data, i, j - 1) * gamm(data, i, j - 1, i, j);
    } else if (j == 0) {
        f = forward(data, i - 1, j) * gamm(data, i - 1, j, i, j);
    } else {
        f = forward(data, i - 1, j) * gamm(data, i - 1, j, i, j) +
            forward(data, i, j - 1) * gamm(data, i, j - 1, i, j) +
            forward(data, i - 1, j - 1) * gamm(data, i - 1, j - 1, i, j);
    }
    // if (i > 1000){
    //         fprintf(stderr, "Got forward value of %12.5e for %d %d \n", f, i, j);
    // }

    F_store[data->N*j + i] = f;
    return f;
}

long double compute_backward(ChannelData *data, int i, int j) {
    // fprintf(stderr, "Computing backwards for %d %d \n", i, j);
    long double b = 0.0;
    if (i == data->N && j == data->Nr) {
        // fprintf(stderr, "Hit Boundary at %d %d\n", i, j );
        b = 1.0;
    } else if (i == data->N) {
        b = 0.0;
    } else if (j == data->Nr) {
        b = backward(data, i + 1, j) * gamm(data, i, j, i + 1, j);
    } else if (i > data->N || j > data->Nr || i < 0 || j < 0) {
        b = 0.0;
    } else {
        b = backward(data, i + 1, j) * gamm(data, i, j, i + 1, j) +
            backward(data, i, j + 1) * gamm(data, i, j, i, j + 1) +
            backward(data, i + 1, j + 1) * gamm(data, i, j, i + 1, j + 1);
    }
    B_store[data->N*j + i] = b;
    return b;
}

long double forward(ChannelData *data, int i, int j) {
    if (F_store[data->N*j + i] != -1.0) {
        return F_store[data->N*j + i];
    }
    if (over_limit(i, j)){
        F_store[data->N*j + i] = 0;
        return 0;
    }
    return compute_forward(data, i, j);
}

long double backward(ChannelData *data, int i, int j) {
    if (B_store[(data->N * j) + i] != -1.0) {
        // fprintf(stderr, "Got from store for %d %d value is %f \n", i, j, B_store[(data->N * j) + i]);
        return B_store[(data->N * j) + i];
    }
    if (over_limit(i, j)){
        B_store[data->N*j + i] = 0;
        // fprintf(stderr, "OL");
        return 0;
    }
    // fprintf(stderr, "Computing");
    return compute_backward(data, i, j);
}

long double likelihood(ChannelData *data, int i, int val){
    // fprintf(stderr, "Attempting likelihood now \n");
    long double l = 0.0;
    for (int j = 0; j < data->Nr; j++){
        if (over_limit(i, j)){
            continue;
        }
        double diag;
        if (data->r[j] == (val + data->w[i]) % 2){
            diag = data->Pt * (1 - data->Ps);
        }
        else{
            diag = data->Pt*data->Ps;
        }

        long double l2 = forward(data, i, j) * diag * backward(data, i + 1, j + 1);
        l += l2;

        long double l3 = forward(data, i, j) * data->Pd * backward(data, i + 1, j);
        l += l3;
    
    }
    l += forward(data, i, data->Nr) * data->Pd * backward(data, i + 1, data->Nr);

    return l;
    
}

int decode(ChannelData *data, FILE* llr_file){
    // fprintf(stderr, "Attempting decode now");
    double *llrs = (double*)malloc(data->Nd * sizeof(double));
     if (llrs==0)
    { fprintf(stderr,"Ran out of memory (while trying to allocate bytes)\n");
        exit(1);
    }
    int successes = 0;
    for (int j = 0; j < data->Nd; j++ ){
        long double pr[2];
        long double p_win = 0.0;
        int ans;
        for (int v = 0; v < 2; v ++){
            pr[v] = likelihood(data, data->info[j], v);
            if (pr[v] > p_win){
                p_win = pr[v];
                ans = v;
            }
            // fprintf(stderr, "Likelihood: %.10f ", pr[v]);
            if (pr[v] == 0){
                printf("Got zero p!");
            }
        }
        if (data->d[j] == ans){
            successes = successes + 1;
        }
        if (pr[0] > 0 && pr[1] > 0){
            llrs[j] = log(pr[0] / pr[1]);
            //printf("LLR %d: %f  \n", j, llrs[j]);
            fprintf(llr_file, "%.6f ", llrs[j]);
            
        }
    }
    fprintf(stderr, "%d/%d successes", successes, data->Nd);
    return 0;
}

int main(int argc, char *argv[]) {
   if (argc != 11) {
        fprintf(stderr, "Usage: %s <priors_file> <llrs_file> <received_file> <info_file> <w_file> <N> <Nr> <Nd> <Pipd> <d_file>\n", argv[0]);
        return 1;
    }
    ChannelData data;
    data.N = atoi(argv[6]);
    data.Nr = atoi(argv[7]);

    F_store = (long double*)malloc((data.N + 1) * (data.Nr + 1) * sizeof(long double));
    B_store = (long double*)malloc((data.N + 1) * (data.Nr + 1) * sizeof(long double));

    initialise_channel_data(&data); 
    char *llr_file = argv[2];
    // printf("Writing llrs to %s", llr_file);
    FILE *llrf;

    llrf = fopen(llr_file, "w");
    if (llrf==NULL)
    { fprintf(stderr,"Can't create file for decoded data: %s\n",llr_file);
        exit(1);
    }
    
    data.Nd = atoi(argv[8]);
    data.Pi = atof(argv[9]);
    data.Pd = atof(argv[9]);
    data.Ps = 0;
    data.Pt = 1.0 - (data.Pi + data.Pd);
    data.priors = d_read_from_file(argv[1], data.Nd);
    data.r = int_read_from_file(argv[3], data.Nr);
    data.w = int_read_from_file(argv[5], data.N);
    data.d = int_read_from_file(argv[10], data.Nd);
    data.info = int_read_from_file(argv[4], data.Nd);

    // for (int i = 0; i < data.Nd; i++){
    //     fprintf(stderr, "Got value %d", data.info[i]);
    // }
    
    // long double h = backward(&data, 5, 5);
    //fprintf(stderr, "Got value %12.5e \n", h);
    
    decode(&data, llrf);

    return 0;

}


