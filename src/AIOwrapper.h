#ifndef AIOWRAPPER_H
#define AIOWRAPPER_H

#include "Definitions.h"
#include "Utility.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <limits>
#include <bitset>
#include <utility>

using namespace std;

typedef struct BufferElement type_buffElement;

struct BufferElement
{
    type_precision* buff;
    int size;
};

typedef struct fileh type_fileh;



struct fileh
{
    string fnameAL;
    string fnameAR;
    string fnameY;

    string io_overhead;


    string fnameOutFiles;
    string fname_dosages;

    int mpi_id;



    list< pair<int,int> >* excl_List;


    bool doublefileType;
    bool fakefiles;
    bool storeBin;

    long double min_p_disp;
    float min_R2_disp;
    bool disp_cov;
    bool storePInd;

    vector< string > ARnames;
    vector< string > Ynames;
    vector< string > ALnames;

    type_precision* Yb;
    type_precision* Ar;
    type_precision* ArDosage;
    type_precision* AL;
    type_precision* B;
    type_buffElement* currentReadBuff;
    type_buffElement* Ar_currentReadBuff;
    list < resultH >* currentWriteBuff;
    int buff_count;

    queue<type_buffElement*> empty_buffers;
    queue<type_buffElement*> full_buffers;

    queue<list < resultH > *> write_empty_buffers;
    queue<list < resultH > *> write_full_buffers;

    queue<type_buffElement*> ar_empty_buffers;
    queue<type_buffElement*> ar_full_buffers;

    //MPI related vars
    int initial_file_pos;

    int index;
    int fileN;
    int fileR;

    int fileM;
    int n;
    int r;
    int l;
    int p;

    int model;
    int modelR;


    int Ar_Amount;
    int Ar_blockSize;
    int Ar_file_blocksize;//used when interactions changes the block size
    int Ar_to_readSize;

    int Y_Amount;
    int y_blockSize;
    int y_to_readSize;

    int b_blockSize;
    int max_b_blockSize;

    bool not_done;
    bool reset_wait;
    bool use_dosages;
    bool add_dosages;
    int dosage_skip;

    bool use_interactions ;
    bool keep_depVar;
    bool use_multiple_interaction_sets;
    string fname_interactions;
    float* interaction_data;
    float* temp_file_data;
    int numInter;
    int ini_IDX_interactions;
    int end_IDX_interactions;

    int seed;
    int Aseed;

    float* dosages;
    vector< vector <float> > cov_2_Terms;
    vector< vector <float> > x_Terms;
    vector< vector <float> > xcov_2_Terms;

    pthread_mutex_t m_more     ;
    pthread_cond_t  condition_more   ;
    pthread_mutex_t m_read     ;

    pthread_cond_t  condition_read  ;
    pthread_mutex_t m_buff_upd   ;
    pthread_mutex_t out_buff_upd   ;

    pthread_t iothread;
    pthread_attr_t attr;

    //barrier
    pthread_barrier_t finalize_barrier;


};

#define fgls_malloc(size) fgls_malloc_impl(__FILE__, __LINE__, size)

enum datatype{ UNSIGNED_SHORT_INT_TYPE = 1,
               SHORT_INT_TYPE,
               UNSIGNED_INT_TYPE,
               INT_TYPE,
               FLOAT_TYPE,
               DOUBLE_TYPE,
               SIGNED_CHAR_TYPE,
               UNSIGNED_CHAR_TYPE };

#define NAMELENGTH 32
#define RESERVEDSPACE 5

typedef struct databel_fvi_header
{
	unsigned short int type;
	unsigned int nelements;
	unsigned int numObservations;
	unsigned int numVariables;
	unsigned int bytesPerRecord;
	unsigned int bitsPerRecord;
	unsigned int namelength;
	unsigned int reserved[RESERVEDSPACE];
} databel_fvi_header;

typedef struct databel_fvi
{
	databel_fvi_header  fvi_header;
	char               *fvi_data;
} databel_fvi;

class AIOwrapper
{
    public:
        AIOwrapper();
        ~AIOwrapper();

        void initialize(struct Settings &params);
        void finalize();


        void load_AL(type_precision** AL);
        void load_ARblock(type_precision** Y, int &blockSize);
        void load_Yblock(type_precision** Y, int &blockSize);
        void reset_Y();
        void reset_AR();

        void getCurrentWriteBuffers(list < resultH >* &sigResults);

        void write_OutFiles(list < resultH >* &sigResults);
        void write_OutSignificant(list < resultH >* sigResults, int min_p_disp, bool disp_cov);

        string io_overhead;
        int realN;



    protected:

    private:

        static void  generate_multinteraction_singleset(float* interaction_data, int cols_data, float* AR, int ar_block_size, int n, float* dest, bool keep_original  );
        static void  generate_singleinteraction_multset(float* interaction_data, int cols_indp_data, float* AR, int ar_block_size, int n, float* dest, bool keep_original  );

        void read_excludeList(list< pair<int,int> >* excl, int &excl_count, int max_excl, string fname_excludeList);
        void read_dosages(string fname_dosages, int expected_count, float* vec);


        void prepare_AR( int desired_blockSize, int n, int totalR, int columnsR);
        void finalize_AR();


        void prepare_Y(int desired_blockSize, int n, int totalY);
        void finalize_Y();

        void prepare_AL( int columns, int n);
        void finalize_AL();

        void prepare_OutFiles();
        void finalize_OutFiles();


        void removeALmissings(list< pair<int,int> >* excl_List,struct Settings params,int &Almissings);
        bool splitpair(int value, list< pair<int,int> >* excl_List,int n);


        static void* async_io(void *ptr );

        struct databel_fvi * load_databel_fvi( const char *path );
        void free_databel_fvi( struct databel_fvi **fvi );
        FILE * fgls_fopen( const char *path, const char *mode );

        void * fgls_malloc_impl( const char* file, long line, size_t size );

        type_fileh FHandler;
        type_fileh* Fhandler;




        databel_fvi* Yfvi;
        databel_fvi* ALfvi;
        databel_fvi* ARfvi;
        databel_fvi* Ifvi;



};

#endif // AIOWRAPPER_H
