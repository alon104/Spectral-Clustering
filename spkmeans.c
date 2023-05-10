# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include "spkmeans.h"

/*Prototypes*/

/*Module Related*/
double** fit_wam(double** mat, int n, int d);
double** fit_ddg(double** mat, int n, int d);
double** fit_gl(double** mat, int n, int d);
double*** fit_jacobi(double** mat, int n);
double** fit_spk(int K, int maxiter, double EPS, int vectorLength,int N, double** clustersArr ,double** dataPointsArr); 

/*Spk related*/
double** wam(double** mat, int n, int d);
double** ddg(double** mat, int n);
double** gl(double** D, double** W, int n);
double*** jacobi(double** L, int n);
double updateOffDiffAndA(double** A, int n, int i, int j, double c, double s);
void updateV(double** V, int n, int i, int j, double c, double s);
int* findPivot(double** A, int n);
double* obtainCS(double** A, int i, int j);
double squaredDistance(double* x, double* y, int n);
double** eye(int n);
void freeMatrix(double** mat, int n);
double** readFile(char* fileName);
void negateColumn(double** matrix, int n, int col);
double off(double** matrix, int n);

/*kmeanspp:*/
struct cord
{
    double value;
    struct cord *next;
};
struct vector
{
    struct vector *next;
    struct cord *cords;
};

int invalidMalloc = 0;
int rows = 0, cols = 0;
int getClosestCentroid(struct vector *curr,struct vector* centroids, int K, int vectorLen);
double distance(struct vector *x, struct vector *y, int vectorLength);
void insertToCluster(struct vector *curr,struct vector* centroids, int index);
void getMeanedCentroid(struct vector *cluster, int *clusterSize, int K, int vectorLen);
int getVectorLength(struct vector vec);
int dataLength(struct vector *vec);
void freeVectors(struct vector* head);
void freeCords(struct cord* head, int vectorLen);
void createNewCentroids(struct vector *newCentroids,int K,int vectorLen);
void zeroClusterSize(int *clusterSize, int K);
int validArg(char* str);
int is_digit(char c);
double** kmeans(int K, int maxiter, double EPS, int vectorLength,int N, double** clustersArr ,double** dataPointsArr); 
void insert_vector(struct vector **head, struct vector **tail);
void insert_cords(double value, struct cord **head, struct cord **tail);
struct vector* fill_dataPoints(double** dataPointsArr,int numOfDataPoints,int vectorLength);
void fillBackCentroids(double** oldCentroids, struct vector *centroidsVectorArr,int K,int vectorLength);
void print_matrix(double** matrix, int n);
void print_jacobi(double** A,double** V,int n);

/*Implementations*/

/*Main*/ 
int main(int argc, char* argv[]) {
    char* goal;
    char* fileName;
    double** res;
    double** matrix;
    int n;
    int d;
    if(argc == 0){
        printf("An Error Has Occurred\n");
        exit(1);
    }
    goal = argv[1];
    fileName = argv[2];
    matrix = readFile(fileName);
    n = rows;
    d = cols;
    if (strcmp(goal, "wam") == 0) {
        res = fit_wam(matrix, n, d);
        print_matrix(res, n);
        freeMatrix(res, n);
        freeMatrix(matrix,n);
        return 0;
    }
    else if (strcmp(goal, "ddg") == 0) {
        res = fit_ddg(matrix, n, d);
        print_matrix(res, n);
        freeMatrix(res, n);
        freeMatrix(matrix,n);
        return 0;
    }
    else if (strcmp(goal, "gl") == 0) {
        res = fit_gl(matrix, n, d);
        print_matrix(res, n);
        freeMatrix(res, n);
        freeMatrix(matrix,n);
        return 0;
    }

    else {  /*Jacobi*/
        double*** AV = fit_jacobi(matrix, n);
        print_jacobi(AV[0], AV[1], n);
        freeMatrix(AV[0], n);
        freeMatrix(AV[1], n);
        free(AV);
        return 0;
        }
    }
    
double** readFile(char* fileName) {
    FILE *fp;
    double **matrix;
    int i, j;
    char ch;

    /*open the file for reading*/
    fp = fopen(fileName, "r");
    if (fp == NULL) {
        printf("An Error Has Occurred\n");
        exit(1);
    }

    /*count the number of rows and columns in the matrix*/
    while ((ch = fgetc(fp)) != EOF) {
        if (ch == ',') {
            cols++;
        }
        else if (ch == '\n') {
            rows++;
            break;
        }
    }
    while ((ch = fgetc(fp)) != EOF) {
        if (ch == '\n') {
            rows++;
        }
    }
    cols++; /*account for the last column*/
    /*allocate memory for the matrix*/
    matrix = (double**)malloc(rows * sizeof(double*));
    if (matrix == NULL) {
        printf("An Error Has Occurred\n");
        exit(1);
    }
    for (i = 0; i < rows; i++) {
        matrix[i] = (double*)malloc(cols * sizeof(double));
        if (matrix[i] == NULL) {
        printf("An Error Has Occurred\n");
        exit(1);
    }
    }
    rewind(fp);
    /*read the matrix elements*/
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            if(fscanf(fp, "%lf,", &matrix[i][j])) {};
        }
        if (fscanf(fp, "\n")) {}; /*read the newline character*/
    }

    /*close the file*/
    fclose(fp);
    return matrix;
}

void print_matrix(double** matrix, int n) {
    int i, j;
    for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
        if (j < n - 1) {
            printf("%.4f,", matrix[i][j]);
        } else {
            printf("%.4f", matrix[i][j]);
        }
    }
    printf("\n");
}

}

void print_jacobi(double** A, double** V, int n){
    int i;
    for (i=0; i < n; i++) {
        if (i < n - 1) {
        printf("%.4f,", A[i][i]);
        } else {
        printf("%.4f", A[i][i]);                    
        }
    }
    printf("\n");
    print_matrix(V, n);
}
/*Exported functions*/
double** fit_wam(double** mat, int n, int d){
    return wam(mat, n, d);
}
double** fit_ddg(double** mat, int n, int d){
    double** W = wam(mat, n, d);
    double** D =  ddg(W, n);
    freeMatrix(W, n);
    return D;
}
double** fit_gl(double** mat, int n, int d){
    double** W = wam(mat, n, d);
    double** D = ddg(W, n);
    double** L = gl(D, W, n);
    freeMatrix(W, n);
    freeMatrix(D, n);
    return L;
}

double*** fit_jacobi(double** mat, int n){
    double*** AV = jacobi(mat, n);
    return AV;
}

double** fit_spk(int K, int maxiter, double EPS, int vectorLength,int N, double** clustersArr ,double** dataPointsArr){
    return kmeans(K, maxiter, EPS, vectorLength, N, clustersArr, dataPointsArr);
}



/*Computing Functions*/
double** wam(double** mat, int n, int d) {
    int i, j;
    double wij;
    double** matrix = (double**)malloc(n * sizeof(double*));
    if (matrix == NULL) {
        printf("An Error Has Occurred\n");
        exit(1);
    }
    for (i = 0; i < n; i++) {
        matrix[i] = (double*)malloc(n * sizeof(double));
        if (matrix[i] == NULL) {
        printf("An Error Has Occurred\n");
        exit(1);
        }
    }
    for (i = 0; i < n; i++) {
        for (j = i; j < n; j++) {
            double dist = squaredDistance(mat[i], mat[j], d);
            wij = i == j ? 0 : exp(-dist/2.0);
            matrix[i][j] = wij;
            matrix[j][i] = wij;
        }
    }
    return matrix;
}

double** ddg(double** matrix, int n) {
    int i;
    int j;
    double di;
    double** diag_matrix = (double**)malloc(n * sizeof(double*));
    if (diag_matrix == NULL) {
        printf("An Error Has Occurred\n");
        exit(1);
    }
    for (i = 0; i < n; i++) {
        diag_matrix[i] = (double*)malloc(n * sizeof(double));
        if (diag_matrix[i] == NULL) {
        printf("An Error Has Occurred\n");
        exit(1);
        }
        di = 0.0;
        for (j = 0; j < n; j++) {
            di += matrix[i][j];
            diag_matrix[i][j] = 0.0;
        }
        diag_matrix[i][i] = di;
    }
    return diag_matrix;
}

double** gl(double** D, double** W, int n) {
    int i;
    int j;
    double** L = (double**)malloc(n * sizeof(double*));
    if (L == NULL) {
        printf("An Error Has Occurred\n");
        exit(1);
    }
    for (i = 0; i < n; i++) {
        L[i] = (double*)malloc(n * sizeof(double));
        if (L[i] == NULL) {
        printf("An Error Has Occurred\n");
        exit(1);
        }
        for (j = 0; j < n; j++) {
            L[i][j] = D[i][j] - W[i][j];
        }
    }
    return L;
}

int* findPivot(double** A, int n) {
    double curr_max_val = -1;
    int curr_max_i = -1;
    int curr_max_j = -1;
    int i, j;
    int* res;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (i == j) {
                continue;
            }
            if (fabs(A[i][j]) > curr_max_val) {
                curr_max_val = fabs(A[i][j]);
                curr_max_i = i;
                curr_max_j = j;
            }
        }
    }
    res = malloc(sizeof(int)*2);
    if (res == NULL) {
        printf("An Error Has Occurred\n");
        exit(1);
    }
    res[0] = curr_max_i;
    res[1] = curr_max_j;
    return res;
}

double* obtainCS(double** A, int i, int j) {
    double theta = (A[j][j] - A[i][i]) / (2 * A[i][j]);
    double sign_theta = (theta >= 0) ? 1 : -1;
    double t = sign_theta / (fabs(theta) + sqrt(pow(theta, 2) + 1));
    double c = 1 / sqrt(pow(t, 2) + 1);
    double s = t * c;
    double* res = malloc(sizeof(double)*4);
    if (res == NULL) {
        printf("An Error Has Occurred\n");
        exit(1);
    }
    res[0] = c;
    res[1] = s;
    res[2] = t;
    res[3] = theta;
    return res;
}

double updateOffDiffAndA(double** A, int n, int i, int j, double c, double s) {
    double off_diff = 0;
    int r;
    double tmp_Aii, tmp_Ajj;
    for (r = 0; r < n; r++){
        if (r != i && r != j) { 
            double tmp_Ari = c * A[r][i] - s * A[r][j];
            double tmp_Arj = c * A[r][j] + s * A[r][i];
            off_diff += 2 * (pow(A[r][i], 2) - pow(tmp_Ari, 2));
            off_diff += 2 * (pow(A[r][j], 2) - pow(tmp_Arj, 2)); /*Symmetric, add to off for both rj and jr. */
            A[r][i] = tmp_Ari;
            A[i][r] = tmp_Ari;
            A[r][j] = tmp_Arj;
            A[j][r] = tmp_Arj;
        }
    }
    tmp_Aii = (pow(c, 2) * A[i][i]) + (pow(s, 2) * A[j][j]) - (2 * s * c * A[i][j]); 
    tmp_Ajj = (pow(s, 2) * A[i][i]) + (pow(c, 2) * A[j][j]) + (2 * s * c * A[i][j]); 
    off_diff += 2 * (pow(A[i][j], 2));
    A[i][i] = tmp_Aii;
    A[j][j] = tmp_Ajj;
    A[i][j] = 0;
    A[j][i] = 0;
    return off_diff;
}

double*** jacobi(double** A, int n) {
    int d;
    double*** result;
    double EPS = 1.0 * pow(10, -5);
    int iter = 0;
    double** V = eye(n); 
    double off_diff = 1;
    while (off_diff > EPS && iter < 100) {
        int* pivot = findPivot(A, n);
        int i = pivot[0], j = pivot[1];
        double* cs = obtainCS(A, i, j);
        double c = cs[0], s = cs[1];
        off_diff = updateOffDiffAndA(A, n, i, j, c, s);
        updateV(V, n, i, j, c, s);
        free(pivot);
        free(cs);
        iter++;
    }
    result = (double***) malloc(sizeof(double**) * 2);
    if (result == NULL) {
        printf("An Error Has Occurred\n");
        exit(1);
    }
    result[0] = A;
    result[1] = V;
    for (d = 0; d < n; d++) {
        if (A[d][d] < 0 && A[d][d] > -0.0001) {
            A[d][d] = 0;
            negateColumn(V, n, d);
        }
    }
    return result;
}

void updateV(double** V, int n, int i, int j, double c, double s) {
    int r;
    double V_ri;
    double V_rj;
    for (r = 0 ; r < n; r++) {
        V_ri = V[r][i] * c - V[r][j] * s;
        V_rj = V[r][i] * s + V[r][j] * c;
        V[r][i] = V_ri;
        V[r][j] = V_rj;
    }
}

double** eye(int n) {
    int i, j;
    double** matrix = (double**) malloc(n * sizeof(double*));
    if (matrix == NULL) {
        printf("An Error Has Occurred\n");
        exit(1);
    }
    for (i = 0; i < n; i++) {
        matrix[i] = (double*) malloc(n * sizeof(double));
        if (matrix[i] == NULL) {
        printf("An Error Has Occurred\n");
        exit(1);
        }
        for (j = 0; j < n; j++) {
            matrix[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }
    return matrix;
}

void freeMatrix(double** mat, int n) {
    int i;
    for (i = 0; i < n; i++) {
        free(mat[i]);
    }
    free(mat);
}


double squaredDistance(double* x, double* y, int n) {
    int i;
    double dist = 0.0;
    for (i = 0; i < n; i++) {
        double diff = x[i] - y[i];
        dist += (diff * diff);
    }
    return dist;
}

void negateColumn(double** matrix, int n, int col) {
    int i;
    for (i = 0; i < n; i++){
        matrix[i][col] *= -1;
    }
}

double off(double** matrix, int n) { 
    int i,j;
    double off = 0;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (j!=i) {
                off += pow(matrix[i][j], 2);
            }
        }
    }
    return off;
}


/*KmeansPP related*/
double** kmeans(int K, int maxiter, double EPS, int vectorLength,int N, double** clustersArr ,double** dataPointsArr){
    struct vector* centroidsVectorArr;
    int i=0;
    int j=0;
    int iter = 0;
    int numOfNewCentroids;
    int *clusterSize;
    struct vector *currVec;
    struct vector *head;
    head = fill_dataPoints(dataPointsArr,N,vectorLength);
    numOfNewCentroids = K;
    clusterSize = malloc(K * sizeof(int));
    centroidsVectorArr = malloc(sizeof(struct vector) * (K));
    if (clusterSize==NULL || centroidsVectorArr==NULL){
       printf("An Error Has Occurred\n");
        exit(1);
    };
    createNewCentroids(centroidsVectorArr, K, vectorLength);
    if (invalidMalloc == 1){
        exit(1);
    };
    for (i=0;i<K;i++){
        struct cord *currCentCords = centroidsVectorArr[i].cords;
            for (j=0; j<vectorLength; j++){
                currCentCords->value = clustersArr[i][j];
                currCentCords = currCentCords->next;}
    } 
    while (iter < maxiter && numOfNewCentroids > 0){
        struct vector* newCentroidsVectorArr = malloc(sizeof(struct vector) * (K));
        struct vector *currVec = head;
        struct vector* tmp;
        if (newCentroidsVectorArr==NULL){
            printf("An Error Has Occurred\n");
            exit(1);
            };
        zeroClusterSize(clusterSize,K);
        createNewCentroids(newCentroidsVectorArr, K , vectorLength);
        numOfNewCentroids = 0;
        while (currVec != NULL){
            int closestIndex = getClosestCentroid(currVec, centroidsVectorArr, K,vectorLength);
            clusterSize[closestIndex]++;
            insertToCluster(currVec, newCentroidsVectorArr, closestIndex);
            currVec = currVec->next;
    }
    getMeanedCentroid(newCentroidsVectorArr,clusterSize, K, vectorLength);
    for (i = 0; i <K; i++){
        if (distance(&newCentroidsVectorArr[i],&centroidsVectorArr[i],vectorLength)>=EPS){
            numOfNewCentroids++;
        }}
        tmp = centroidsVectorArr;
        centroidsVectorArr = newCentroidsVectorArr;
        newCentroidsVectorArr = tmp;
        for (i=0; i<K;i++){
            freeCords(newCentroidsVectorArr[i].cords,vectorLength);
        }
    
    free(newCentroidsVectorArr);
    iter++;
}

fillBackCentroids(clustersArr,centroidsVectorArr, K,vectorLength);

for (i=0; i<K;i++){
    freeCords(centroidsVectorArr[i].cords,vectorLength);
    }
currVec = head;
free(centroidsVectorArr);
free(clusterSize);
for (i = 0; i <N; i++){
    freeCords(currVec->cords,vectorLength);
    currVec = currVec->next;
}
freeVectors(head->next);
free(head);
return clustersArr;
}

int getClosestCentroid(struct vector *curr,struct vector* centroids, int K, int vectorLen){
    int closestIndex=0;
    double minimum = distance(curr,&centroids[0],vectorLen);
    int i;
    for (i=0;i<K;i++){
        double dist = distance(curr,&centroids[i],vectorLen);
        if (dist < minimum){
            minimum = dist;
            closestIndex = i;
        }
    }
    return closestIndex;
}

double distance(struct vector *x, struct vector *y, int vectorLength){
    struct cord *yCord = y->cords;
    struct cord *xCord = x->cords;
    double sum = 0.0;
    int i=0;
    for (i = 0; i < vectorLength;i++){
        sum += pow((yCord->value-xCord->value),2);
        yCord = yCord->next;
        xCord = xCord->next;
    }
    return sqrt(sum);
}

void insertToCluster(struct vector *curr,struct vector* centroids, int index){
    struct cord *tmp; 
    struct cord *centroidCords;
    tmp = curr->cords;
    centroidCords = centroids[index].cords;
    while (tmp != NULL){
        centroidCords->value+=tmp->value;
        centroidCords = centroidCords->next;
        tmp = tmp->next;
    }
}

void getMeanedCentroid(struct vector *cluster, int *clusterSize, int K, int vectorLen){
    int i;
    int j;
    for (i = 0; i <K; i++){
        struct cord *curr = cluster[i].cords; 
        for (j = 0; j < vectorLen; j++){
            curr->value/=clusterSize[i];
            curr = curr->next;
        }
    }
}

int getVectorLength(struct vector vec){
    struct cord *pt;
    int len = 0;
    pt = vec.cords;
    while (pt!=NULL){
        len++;
        pt = pt->next;
    }
    return len;
}

int dataLength(struct vector *vec){
    struct vector *curr;
    int len = 0;
    curr = vec;
    while (curr!=NULL){
        len++;
        curr = curr->next;
    }
    return len;
}

void freeVectors(struct vector* headVec)
{
   struct vector* tmp;

   while (headVec != NULL)
    {
               
       tmp = headVec;
       headVec = headVec->next;
       free(tmp);

       
    }

}
void freeCords(struct cord* headCord, int vectorLen)
{
    struct cord* tmp;
    int i;
   for (i = 0; i < vectorLen; i++)
    {
       tmp = headCord;
       headCord = headCord->next;
       free(tmp);
    }
    free(headCord);
}

void createNewCentroids(struct vector *newCentroids,int K,int vectorLen){
    int i;
    int j;
    struct cord *curr;
    for (i=0;i<K;i++){
        newCentroids[i].cords = malloc(sizeof(struct cord));
        if (newCentroids[i].cords==NULL){
            printf("An Error Has Occurred\n");
            invalidMalloc = 1;
            return;
            };
        curr = newCentroids[i].cords;
        for (j=0;j<vectorLen;j++){
            curr->value = 0;
            curr->next = malloc(sizeof(struct cord));
            if (curr->next==NULL){
            printf("An Error Has Occurred\n");
            invalidMalloc = 1;
            return;
            };
            curr = curr->next;
        }
    }
    } 

void zeroClusterSize(int *clusterSize, int K){
    int i;
    for (i = 0; i < K; i++){
        clusterSize[i] = 0;
    }
}

int validArg(char* str){
    int i=0;
    char curr = *(str);
    while (curr){
        if (!is_digit(curr)){
            return 0;
        }
        i++;
        curr = *(str + i);
    }
    return 1;
}

int is_digit(char c){
    if (c >= '0' && c <= '9'){
        return 1;
    }
    else{
        return 0;
    }
}

void insert_vector(struct vector **head, struct vector **tail) {
  struct vector *new_vec = malloc(sizeof(struct vector));
  if(new_vec == NULL){
    printf("An Error Has Occurred\n");
    exit(1);
  }
  new_vec->cords = NULL;
  new_vec->next = NULL;

  if (*head == NULL) {
    *head = new_vec;
  } else {
    (*tail)->next = new_vec;
  }
  *tail = new_vec;
}

void insert_cords(double value, struct cord **head, struct cord **tail) {
  struct cord *new_cords = malloc(sizeof(struct cord));
  if(new_cords == NULL){
    printf("An Error Has Occurred\n");
    exit(1);
  }
  new_cords->value = value;
  new_cords->next = NULL;

  if (*head == NULL) {
    *head = new_cords;
  } else {
    (*tail)->next = new_cords;
  }
  *tail = new_cords;
}

struct vector* fill_dataPoints(double** dataPointsArr,int numOfDataPoints,int vectorLength){
    int i, j;
    struct cord *cords_tail;
    struct vector *vec_head = NULL, *vec_tail = NULL;
    for (i = 0; i < numOfDataPoints; i++) {
        insert_vector(&vec_head, &vec_tail);
        cords_tail = NULL;
        for (j = 0; j < vectorLength; j++) {
            insert_cords(dataPointsArr[i][j], &(vec_tail->cords), &cords_tail);
  }
}
  return vec_head;
}

void fillBackCentroids(double** oldCentroids, struct vector *centroidsVectorArr,int K,int vectorLength){
    int i=0;
    int j=0;
    struct cord *curr_cord;
    for (i = 0; i < K; i++){
        curr_cord = centroidsVectorArr[i].cords;
        for (j = 0; j < vectorLength; j++) {
            oldCentroids[i][j] = curr_cord->value;
            curr_cord = curr_cord->next;
        }
    }
}

