/* ---------------------------------------------------------------
Práctica 1.
Código fuente : nbfast.c
Grau Informàtica
Rafel Salgueiro Santacreu 48257208B
--------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <unistd.h>
#include <pthread.h>
#include <string.h>
#include <semaphore.h>

#ifdef D_GLFW_SUPPORT
    #include<GLFW/glfw3.h>
#endif

// Macros to make code a little bit easier to understand because for speedup reasons, I'll use only 1D arrays
#define PX(i) (3*i+1)
#define PY(i) (3*i+2)
#define MASS(i) (3*i+3)

#define VX(i) (4*i+0)
#define VY(i) (4*i+1)
#define AX(i) (4*i+2)
#define AY(i) (4*i+3)

double G=0.0001;
double dt=0.005;
double rcutoff=0.35;
double rlimit=0.03;

struct Node{
    struct Node *children[4];
    int external;

    double CMX;
    double CMY;
    double mass;
    double TRX;
    double TRY;

    double LLX;
    double LLY;

    double GCX;
    double GCY;
};

struct globalVariables{
    struct Node* tree;
    int nShared;
	int steps;
    double *sharedBuff;
    double *localBuff;
    int *indexes;
    int possibleThreads;
    int id;
};
typedef struct globalVariables globalVariables, *globalVariablesPtr;

struct threadBT{
    struct Node* node;
    double* shrdBuff;
    int *indexes;
    int n;
    int Threads;
};
typedef struct threadBT threadBT, *threadBTptr;



pthread_t *threads;
globalVariablesPtr globalVars;

//Variables de sincronización
sem_t initCalculateForce;
sem_t endCalculateForceAndMoveParticle;
pthread_barrier_t itBarrier;
pthread_barrier_t CalculateForceEndBarrier;
pthread_mutex_t staticsMutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t globalVariablesMutex = PTHREAD_MUTEX_INITIALIZER;




void *buildTreeThread (threadBTptr job);

void buildTree(struct Node* node, double* shrdBuff, int *indexes, int n, int MThreads){
    pthread_t *MiTids = NULL;
    threadBTptr MiJobs = NULL;
    int Jobs, NewThreads = 0, RemainingThreads = 0, PendingThreads = 0, CurrentJob = 0;

    if(n==1){ //This is an external node!
        node->external=1;
        node->CMX=shrdBuff[PX(indexes[0])];
        node->CMY=shrdBuff[PY(indexes[0])];
        node->mass=shrdBuff[MASS(indexes[0])];
    } else {
        node->external=0;
		//Arrays of indexes of particles per quartile
        int *NEi = (int *) malloc(sizeof(int)*n);
        int *NWi = (int *) malloc(sizeof(int)*n);
        int *SWi = (int *) malloc(sizeof(int)*n);
        int *SEi = (int *) malloc(sizeof(int)*n);
        int NWc=0, SWc=0,SEc=0, NEc=0;
        int i;
 
		/** For each particle we will check where is it located relative to the geometric center,
			to sort them into the 4 children nodes**/
        for(i=0;i<n;i++){
            if(shrdBuff[PY(indexes[i])] < node->GCY ){ //South half
                if(shrdBuff[PX(indexes[i])] < node->GCX){ //West wing
                    SWi[SWc]=indexes[i];
                    SWc++;
                } else {
                    SEi[SEc]=indexes[i];
                    SEc++;
                }
            } else { //North half
                if(shrdBuff[PX(indexes[i])] < node->GCX){ //West wing
                    NWi[NWc]=indexes[i];
                    NWc++;
                } else {
                    NEi[NEc]=indexes[i];
                    NEc++;
                }
            }
        }
        if (MThreads > 0){              //si hay hilos disponibles

            Jobs=0;
            //Contamos trabajos a realizar
            if (NEc > 0) Jobs++;
            if (NWi > 0) Jobs++;
            if (SWc > 0) Jobs++;
            if (SEc > 0) Jobs++;

            if (MThreads > Jobs){       //si hay mas hilos que trabajos, se crean tantos hilos como trabajos
                NewThreads = Jobs;
            } else {                    //si hay mas trabajos que hilos, se crean tantos hilos como se puedan
                NewThreads = MThreads;
            }
            RemainingThreads = MThreads - Jobs;         //se calculan los hilos sobrantes
            if (RemainingThreads < 0) RemainingThreads = 0;

            MiTids = malloc (sizeof(pthread_t)*NewThreads);
            if (MiTids == NULL) {
                perror("Error reservar Tids vector.");
            }   
            MiJobs = malloc (sizeof(threadBT)*NewThreads);
            if (MiJobs == NULL) {
                perror("Error reservar Jobs vector.");
            }
            PendingThreads = NewThreads;
            CurrentJob = 0;            
        }
		//If there are particles in the NorthWest quarter
        if(NEc>0){
			//This instruction declares a new node in the position 0
            node->children[0]= malloc(sizeof *node->children[0]);
			//We give the values of the Low Left and Top Right corner, and also the geometric center.
            node->children[0]->TRX=node->TRX;
            node->children[0]->TRY=node->TRY;
            node->children[0]->LLX=node->GCX;
            node->children[0]->LLY=node->GCY;
            node->children[0]->GCX=(node->GCX+node->TRX)/2;
            node->children[0]->GCY=(node->GCY+node->TRY)/2;
            if (PendingThreads > 0){                                    //si quedan hilos disponibles
                MiJobs[CurrentJob].node = node->children[0];            //se le asigna un trabajo
                MiJobs[CurrentJob].shrdBuff = shrdBuff;
                MiJobs[CurrentJob].indexes = NEi;
                MiJobs[CurrentJob].n = NEc;
                MiJobs[CurrentJob].Threads = RemainingThreads/PendingThreads;       //se le asignan los hilos sobrantes
                if(pthread_create(&(MiTids[CurrentJob]), NULL, (void *(*)(void *)) buildTreeThread, (void*) &(MiJobs[CurrentJob]))!=0)
                    perror("Error crear hilo");
                RemainingThreads -= MiJobs[CurrentJob].Threads;
                PendingThreads--;
                CurrentJob++;
            } else {                                                    //si no quedan hilos disponibles
                buildTree(node->children[0], shrdBuff, NEi, NEc, PendingThreads);
            }
        } else {
			//If not, we set the children to null
            node->children[0]=NULL;
        }
                                          //se desbloquea el mutex
		//The next three blocks are exactly the same thing but for the other three nodes
        if(NWc>0){
            node->children[1]= malloc(sizeof *node->children[1]);
            node->children[1]->TRX=node->GCX;
            node->children[1]->TRY=node->TRY;
            node->children[1]->LLX=node->LLX;
            node->children[1]->LLY=node->GCY;
            node->children[1]->GCX=(node->LLX+node->GCX)/2;
            node->children[1]->GCY=(node->GCY+node->TRY)/2;
            if (PendingThreads > 0){                                //si hay hilos disponibles
                MiJobs[CurrentJob].node = node->children[1];        //se asigna trabajo a hilo
                MiJobs[CurrentJob].shrdBuff = shrdBuff;
                MiJobs[CurrentJob].indexes = NWi;
                MiJobs[CurrentJob].n = NWc;
                MiJobs[CurrentJob].Threads = RemainingThreads/PendingThreads;
                if(pthread_create(&(MiTids[CurrentJob]), NULL, (void *(*)(void *)) buildTreeThread, (void*) &(MiJobs[CurrentJob]))!=0)
                    perror("Error crear hilo");

                RemainingThreads -= MiJobs[CurrentJob].Threads;
                PendingThreads--;
                CurrentJob++;
            } else {
                buildTree(node->children[1], shrdBuff, NWi, NWc, PendingThreads);
            }
        } else {
            node->children[1]=NULL;
        }
        if(SWc>0){
            node->children[2]= malloc(sizeof *node->children[2]);
            node->children[2]->TRX=node->GCX;
            node->children[2]->TRY=node->GCY;
            node->children[2]->LLX=node->LLX;
            node->children[2]->LLY=node->LLY;
            node->children[2]->GCX=(node->LLX+node->GCX)/2;
            node->children[2]->GCY=(node->LLY+node->GCY)/2;
            if (PendingThreads > 0){                                //si hay hilos disponibles
                MiJobs[CurrentJob].node = node->children[2];        //se le asigna un trabajo
                MiJobs[CurrentJob].shrdBuff = shrdBuff;
                MiJobs[CurrentJob].indexes = SWi;
                MiJobs[CurrentJob].n = SWc;
                MiJobs[CurrentJob].Threads = RemainingThreads/PendingThreads;
                if(pthread_create(&(MiTids[CurrentJob]), NULL, (void *(*)(void *)) buildTreeThread, (void*) &(MiJobs[CurrentJob]))!=0)
                    perror("Error crear hilo");
                RemainingThreads -= MiJobs[CurrentJob].Threads;
                PendingThreads--;
                CurrentJob++;
            } else {
                buildTree(node->children[2], shrdBuff, SWi, SWc, PendingThreads);
            }
        } else {
            node->children[2]=NULL;
        }
        if(SEc>0){
            node->children[3]= malloc(sizeof *node->children[3]);
            node->children[3]->TRX=node->TRX;
            node->children[3]->TRY=node->GCY;
            node->children[3]->LLX=node->GCX;
            node->children[3]->LLY=node->LLY;
            node->children[3]->GCX=(node->GCX+node->TRX)/2;
            node->children[3]->GCY=(node->LLY+node->GCY)/2;
            if (PendingThreads > 0){                                    //Si hay hilos disponibles
                MiJobs[CurrentJob].node = node->children[3];            //Se le asigna el trabajo
                MiJobs[CurrentJob].shrdBuff = shrdBuff;
                MiJobs[CurrentJob].indexes = SEi;
                MiJobs[CurrentJob].n = SEc;
                MiJobs[CurrentJob].Threads = RemainingThreads/PendingThreads;
                if(pthread_create(&(MiTids[CurrentJob]), NULL, (void *(*)(void *)) buildTreeThread, (void*) &(MiJobs[CurrentJob]))!=0)
                    perror("Error crear hilo");
                RemainingThreads -= MiJobs[CurrentJob].Threads;
                PendingThreads--;
                CurrentJob++;
            } else {                                                       //Si no hay hilos disponibles
                buildTree(node->children[3], shrdBuff, SEi, SEc, PendingThreads);
            }
        } else {
            node->children[3]=NULL;
        }
        node->mass=0;
        node->CMX=0;
        node->CMY=0;
		//Now that we have finished building the 4 trees beneath this node, we calculate the Center of Mass
		//based on the center of mass of the children
        for(i=0;i<4;i++){
            if(node->children[i]!=NULL){
                node->mass+=node->children[i]->mass;
                node->CMX+=node->children[i]->CMX*node->children[i]->mass;
                node->CMY+=node->children[i]->CMY*node->children[i]->mass;
            }
        }
        node->CMX=node->CMX/node->mass;
        node->CMY=node->CMY/node->mass;
		//And tadaaa
    }
    for (int c = 0; c < NewThreads; c++) {
        if(pthread_join(MiTids[c], NULL)!=0) perror("Error join");      //Esperamos a que terminen los hilos
    }
}

void *buildTreeThread (threadBTptr job){
    buildTree(job->node, job->shrdBuff, job->indexes, job->n, job->Threads);
    return(NULL);
}

void calculateForce(struct Node *tree, double *shrdBuff, double *localBuff, int index){
    double distance = sqrt((tree->CMX-shrdBuff[PX(index)])*(tree->CMX-shrdBuff[PX(index)])+
                           (tree->CMY-shrdBuff[PY(index)])*(tree->CMY-shrdBuff[PY(index)]));
	//First we check if the node is not actually the same particle we are calculating
    if(distance>0){
		//Now, we know it is not because the is some distance between the Center of Mass and the particle
		//If the node is external (only contains one particle) or is far away enough, we calculate the force with the center of mass
        if(distance>rcutoff || tree->external){
            double f;
            if(distance<rlimit){
                f=G*tree->mass/(rlimit*rlimit*distance);
            } else {
                f=G*tree->mass/(distance*distance*distance);
            }
            localBuff[AX(index)]+=f*(tree->CMX-shrdBuff[PX(index)]);
            localBuff[AY(index)]+=f*(tree->CMY-shrdBuff[PY(index)]);
        } else {
			//If not, we recursively call the calculateForce() function in the children that are not empty.
            int i;
            for(i=0;i<4;i++){
                if(tree->children[i]!=NULL){
                    calculateForce(tree->children[i],shrdBuff,localBuff,index);
                }
            }
        }
    }
}

void moveParticle(double *shrdBuff, double *localBuff, int index){
    //Unprecise but fast euler method for solving the time differential equation
	double oldX=shrdBuff[PX(index)];
    double oldY=shrdBuff[PY(index)];
    shrdBuff[PX(index)]+=localBuff[VX(index)]*dt+localBuff[AX(index)]*dt*dt*0.5;
    shrdBuff[PY(index)]+=localBuff[VY(index)]*dt+localBuff[AY(index)]*dt*dt*0.5;
    localBuff[VX(index)]=(shrdBuff[PX(index)]-oldX)/dt;
    localBuff[VY(index)]=(shrdBuff[PY(index)]-oldY)/dt;
}

#ifdef D_GLFW_SUPPORT
void drawParticle(double *shrdBuff, double *radius, int index){
    glBegin(GL_TRIANGLE_FAN);
    int k;
    glVertex2f(shrdBuff[PX(index)],shrdBuff[PY(index)]);
    for(k=0;k<20;k++){
        float angle=(float) (k)/19*2*3.141592;
        glVertex2f(shrdBuff[PX(index)]+radius[index]*cos(angle),shrdBuff[PY(index)]+radius[index]*sin(angle));
    }
    glEnd();
}

void drawBarnesHutDivisions(struct Node *rootNode){
    if(!rootNode->external){
        glBegin(GL_LINES);
        glVertex2f(rootNode->GCX,rootNode->LLY);
        glVertex2f(rootNode->GCX,rootNode->TRY);
        glVertex2f(rootNode->LLX,rootNode->GCY);
        glVertex2f(rootNode->TRX,rootNode->GCY);
        glEnd();
        int i;
        for(i=0;i<4;i++){
            if(rootNode->children[i]!=NULL){
                drawBarnesHutDivisions(rootNode->children[i]);
            }
        }
    }
}
#endif

void SaveGalaxy(int count, int nShared, int *indexes, double *sharedBuff);
void SaveGalaxyFile(char *filename, int nShared, int *indexes, double *sharedBuff);

void SaveGalaxy(int count, int nShared, int *indexes, double *sharedBuff)
{
    char filename[100];
    sprintf(filename,"./res/galaxy_%dB_%di.out",nShared,count);
    SaveGalaxyFile(filename, nShared, indexes, sharedBuff);
}

void SaveGalaxyFile(char *filename, int nShared, int *indexes, double *sharedBuff)
{
    int i;
    FILE *res = fopen(filename,"w");

    fprintf(res,"%d\n",nShared);
    for(i=0;i<nShared;i++){
        fprintf(res,"%d\t%e\t%e\t%e\n",indexes[i],sharedBuff[PX(indexes[i])],sharedBuff[PY(indexes[i])],sharedBuff[MASS(indexes[i])]);
    }
    fclose(res);
}

void ReadGalaxyFile(char *filename, int *nShared, int **indexes, double **sharedBuff)
{
    int i,ind;
    FILE *input;

    printf("Reading file %s\n",filename);

    input = fopen(filename,"r");
    if (input==NULL) {
        printf("Error opening file.\n");
        exit(1);
    }
    // Read number of bodies.
    if (fscanf(input,"%d\n",nShared)<1){
        printf("Error reading number of particles.\n");
        exit(1);
    }


    printf("Reading %d bodies\n",*nShared);

    // Reserve memory for indexes and particles.
    *indexes = (int*) malloc(sizeof(int)*(*nShared));
    *sharedBuff = (double *) malloc(sizeof(double)*(3*(*nShared)+1));

    for(i=0;i<(*nShared);i++){
        (*indexes)[i]=i;
    }

    for(i=0;i<(*nShared);i++){
        if (fscanf(input,"%d\t%le\t%le\t%le\n", &ind,&((*sharedBuff)[PX((*indexes)[i])]),&((*sharedBuff)[PY((*indexes)[i])]),&((*sharedBuff)[MASS((*indexes)[i])]))<4){
            printf("Error reading number of particles.\n");
            exit(1);
        }
        //printf("Body %d: (%le,%le) %le\n", ind,(*sharedBuff)[PX((*indexes)[i])],(*sharedBuff)[PY((*indexes)[i])],(*sharedBuff)[MASS((*indexes)[i])]);
    }

    fclose(input);
}


void threadFunction(globalVariablesPtr gV){
    //unpack global variables
    pthread_mutex_lock(&globalVariablesMutex);
    struct Node* tree;
    int id = gV->id;
    int possibleThreads = gV->possibleThreads;
    int nShared = gV->nShared;
    int steps = gV->steps;
    int *indexes;
    double *sharedBuff;
    double *localBuff;

    int countIteration = 0;
    int from;
    int to;
    int particlesPerThread = nShared/possibleThreads;
    from = (id-1) * particlesPerThread;
    to = from + particlesPerThread;
    if (id == possibleThreads){
        to = nShared;
    }
    pthread_mutex_unlock(&globalVariablesMutex);

    pthread_barrier_wait(&itBarrier);               //Wait all the threads to go together

    while (1){  
        sem_wait(&initCalculateForce);                      //Wait for the main thread to finish initializing the tree
        //unpack global variables
        pthread_mutex_lock(&globalVariablesMutex);
        localBuff = gV->localBuff;
        indexes = gV->indexes;
        sharedBuff = gV->sharedBuff;
        tree = gV->tree;
        pthread_mutex_unlock(&globalVariablesMutex);

        int i = 0;
        for(i=from; i < to; i++){
            //Set initial accelerations to zero
            localBuff[AX(indexes[i])]=0;
            localBuff[AY(indexes[i])]=0;

            int s;
            for(s=0;s<4;s++){
                //Recursively calculate accelerations
                if(tree->children[s]!=NULL){
                    // If there are free threads to be created, we execute the next recursive call concurrently.

                    calculateForce(tree->children[s], sharedBuff, localBuff, indexes[i]);

                }
            }

            //Calculate new position
            moveParticle(sharedBuff,localBuff,indexes[i]);

            if (sharedBuff[PX(indexes[i])]<=0 || sharedBuff[PX(indexes[i])]>=1 || sharedBuff[PY(indexes[i])] <=0 || sharedBuff[PY(indexes[i])] >= 1) {
                // If the particle is out of the limits, we count it for the statistics.
                // removedParticles++;
            }
        }

        pthread_mutex_lock(&globalVariablesMutex);
        countIteration++;
        pthread_mutex_unlock(&globalVariablesMutex);

        int blockedThreads = pthread_barrier_wait(&CalculateForceEndBarrier);
        if (blockedThreads == PTHREAD_BARRIER_SERIAL_THREAD){
            sem_post(&endCalculateForceAndMoveParticle);
        } 
        
        if (countIteration >= steps){
            pthread_exit(NULL);
        }
    }
}

#define DSaveIntermediateState 1
#define DIntervalIntermediateState 100
#define DShowStatistics 1
#define DIntervalStatistics 1

clock_t StartTime, EndTime;
double TimeSpent;

void ShowWritePartialResults(int count,int nOriginal, int nShared, int *indexes, double *sharedBuff)
{
    if (DSaveIntermediateState && !(count % DIntervalIntermediateState))
        SaveGalaxy(count, nOriginal, indexes, sharedBuff);

    if (DShowStatistics && !(count % DIntervalStatistics))
    {
        int i=0;
        double CurrentTime;
        CurrentTime = clock();
        TimeSpent = (double)(CurrentTime - StartTime) / CLOCKS_PER_SEC;
        //Mins = (int)TimeSpent/60;
        //Secs = (TimeSpent-(Mins*60));
        printf("[%.3f] Iteration %d => %d Bodies (%d) \t(Body %d: (%le, %le) %le).\n",TimeSpent, count, nShared, nOriginal, i, sharedBuff[PX(indexes[i])],sharedBuff[PY(indexes[i])],sharedBuff[MASS(indexes[i])]);
    }
}



int main(int argc, char *argv[]){
    int nShared=500;
	int steps=100;
    double *sharedBuff;
    double *localBuff;
    double *radius;
    int *indexes, i;
    char filename[100];
    int possiblePthreads = 4;
    
    //inicializacion de las variables de condicion
    sem_init (&initCalculateForce, 0, 0);
    sem_init (&endCalculateForceAndMoveParticle, 0, 0);
    pthread_barrier_init(&itBarrier, NULL, possiblePthreads);
    pthread_barrier_init(&CalculateForceEndBarrier, NULL, possiblePthreads);



    printf("NBody with %d arguments.\n",argc);
    StartTime = clock();


	if(argc>1){
		nShared=atoi(argv[1]);
		if(argc>2){
			steps=atoi(argv[2]);
		}
	}

    if(argc>3 && access(argv[3], F_OK) == 0)
    {
        printf("Read file..\n");
        /* Read bodies initial state from file */
        ReadGalaxyFile(argv[3], &nShared, &indexes, &sharedBuff);
        argc--;
    }
    else
    {   /* Inicialize the bodies randomly */

        //Buffers to hold the position of the particles and their mass
        sharedBuff = (double *) malloc(sizeof(double) * (3 * nShared + 1));

        srand(time(NULL));
        for (i = 0; i < nShared; i++) {
            //I start with an almost random distribution of particles
            sharedBuff[PX(i)] = (float) (i) / (nShared - 1) * 0.8 + 0.1;
            sharedBuff[PY(i)] = (float) (rand() % 4096) / 4095 * 0.8 + 0.1;
            //With a random Mass between 1 and 3
            sharedBuff[MASS(i)]=(double) (rand()%2048)/2047*2+1;
        }

        //Index array, to speed up the creation of the tree (faster than passing the 3 floats per particle of x,y and mass)
        indexes = (int*) malloc(sizeof(int)*nShared);
        for(i=0;i<nShared;i++){
            indexes[i]=i;
        }
    }

    if (argc>4)         //Si se le assigna un nombre de threads
    {
        possiblePthreads = atoi(argv[4]);
    }
    printf("NBody with %d threads.\n",possiblePthreads);

    int nLocal=nShared;
    int nOriginal=nShared;
    //Buffer to hold velocity in x and y, and acceleration in x and y also
    localBuff = (double *) malloc(sizeof(double)*(4*nLocal));
    //This is for opengl
    radius = (double *) malloc(sizeof(double)*(nShared));

    for(i=0;i<nShared;i++){
        // init bodies mass
        radius[i]=sqrt(sharedBuff[MASS(i)])*0.0025;

        //With zero speed, and zero acceleration
        localBuff[VX(i)]=0;
        localBuff[VY(i)]=0;
        localBuff[AX(i)]=0;
        localBuff[AY(i)]=0;
    }

    //This is the main node, the one that holds the first four children nodes that make the calculation zone
    struct Node* tree = malloc(sizeof *tree);
	//LLX is the x coordinate of the Low Left corner
    tree->LLX=0;
	//This is the y coordinate..
    tree->LLY=0;

	//Now the same but for the top right corner
    tree->TRX=1;
    tree->TRY=1;
	//The coordinates of the geometric center of the node in x and y
    tree->GCX=0.5;
    tree->GCY=0.5;

    // Save initial state.
    sprintf(filename,"./res/galaxy_%dB_initial.out",nOriginal);
    SaveGalaxyFile(filename, nShared, indexes, sharedBuff);

    threads = malloc(sizeof(pthread_t) * possiblePthreads);
    if (threads == NULL){
        perror("Error allocating memory for threads");
        exit(1);
    }
    globalVars = malloc(sizeof(globalVariables) * possiblePthreads);
    if (globalVars == NULL){
        perror("Error allocating memory for global variables");
        exit(1);
    }


    for (int i = 0; i < possiblePthreads; i++){
        globalVars[i].possibleThreads = possiblePthreads;
        globalVars[i].nShared = nShared;
        globalVars[i].localBuff = localBuff;
        globalVars[i].indexes = indexes;
        globalVars[i].sharedBuff = sharedBuff;
        globalVars[i].tree = tree;
        globalVars[i].steps = steps;
        globalVars[i].id = i+1;

        if (pthread_create(&threads[i], NULL, (void *(*)(void *)) threadFunction, (void *)&globalVars[i])){
            perror("Error creating threads");
        }
    }
    int count=1;
	//If we need to visualize
#ifdef D_GLFW_SUPPORT
	if(argc>4){
		//If you only care about the algorithm, skip until next comment
	    if(!glfwInit()){
    	    printf("Failed to start GLFW\n");
        	return -1;
    	}
    	GLFWwindow *window = glfwCreateWindow(2000,2000,"Simulation",NULL,NULL);
    	if(!window){
        	printf("Failed to open window\n");
        	return -1;
    	}
    	glfwMakeContextCurrent(window);
    	glfwSwapInterval(1);

    	glMatrixMode(GL_PROJECTION);
    	glLoadIdentity();
    	glOrtho(0,1,0,1,0,1);
    	glMatrixMode(GL_MODELVIEW);

    	while(!glfwWindowShouldClose(window) && count<=steps){
        	glClear(GL_COLOR_BUFFER_BIT);

			double t=glfwGetTime();
			//We build the tree, which needs a pointer to the initial node, the buffer holding position and mass of the particles, indexes and number of particles
        	buildTree(tree,sharedBuff,indexes,nShared, possiblePthreads);
        	//Now that it is built, we calculate the forces per particle
            for (int i = 0; i < possiblePthreads; i++){
                globalVars[i].localBuff = localBuff;
                globalVars[i].indexes = indexes;
                globalVars[i].sharedBuff = sharedBuff;
                globalVars[i].tree = tree;
            }
            for (int a = 0; a < possiblePthreads; a++){
                sem_post(&initCalculateForce);
            }

            sem_wait(&endCalculateForceAndMoveParticle);

            tree = globalVars[0].tree;
            localBuff = globalVars[0].localBuff;
            indexes = globalVars[0].indexes;
            sharedBuff = globalVars[0].sharedBuff;

            SaveGalaxy(count, nShared, indexes, sharedBuff);

			//This is only for visualization
        	drawBarnesHutDivisions(tree);
        	int k;
        	for(k=0;k<nShared;k++){
            	drawParticle(sharedBuff,radius,indexes[k]);
        	}

			t=glfwGetTime()-t;
			if(t<0.013){
				usleep(1000*1000*(0.013-t));
			}

        	glfwSwapBuffers(window);
        	glfwPollEvents();

            ShowWritePartialResults(count, nOriginal, nShared, indexes, sharedBuff);

            //We advance one step
			count++;
    	}
    	glfwTerminate();
	} else {
#endif
		//This is the pure algorithm, without visualization
		//system("mkdir res");
    	while(count<=steps){
			//First we build the tree
        	buildTree(tree,sharedBuff,indexes,nShared,possiblePthreads);

            for (int i = 0; i < possiblePthreads; i++){
                globalVars[i].localBuff = localBuff;
                globalVars[i].indexes = indexes;
                globalVars[i].sharedBuff = sharedBuff;
                globalVars[i].tree = tree;
            }
            for (int a = 0; a < possiblePthreads; a++){
                sem_post(&initCalculateForce);
            }

            sem_wait(&endCalculateForceAndMoveParticle);

            tree = globalVars[0].tree;
            localBuff = globalVars[0].localBuff;
            indexes = globalVars[0].indexes;
            sharedBuff = globalVars[0].sharedBuff;
			//To be able to store the positions of the particles
            ShowWritePartialResults(count,nOriginal, nShared, indexes, sharedBuff);
            //We advance one step
			count++;
		}
#ifdef D_GLFW_SUPPORT
	}
#endif

    EndTime = clock();
    TimeSpent = (double)(EndTime - StartTime) / CLOCKS_PER_SEC;
    printf("NBody Simulation took %.3f seconds.\n",TimeSpent);

    // Save initial state.
    sprintf(filename,"./res/galaxy_%dB_%di_final.out",nOriginal, count-1);
    SaveGalaxyFile(filename, nShared, indexes, sharedBuff);

    // Free memory.
    sem_destroy(&initCalculateForce);
    sem_destroy(&endCalculateForceAndMoveParticle);
    pthread_barrier_destroy(&itBarrier);
    pthread_barrier_destroy(&CalculateForceEndBarrier);
    pthread_mutex_destroy(&staticsMutex);
    pthread_mutex_destroy(&globalVariablesMutex);


	free(sharedBuff);
	free(localBuff);
	free(radius);
	free(indexes);

    return 0;
}