#ifndef _FIRST_
#define _FIRST_
struct pHDU {
    char *observer, *project, *telescope;
    double ant_x,ant_y,ant_z;
    char *fend;
    int ibeam;
    int numChan, fd_had;
    char *fd_pol;
    float fd_sang, fd_xy;
    char *bend, *befcfg;
    int bep, bedcc, obsfch;
    float bede,obscf,obsbw,dm;
    double tcycle;
    char *obsm, *dobs;
    char *pnid, *reffrm;
    char srcn[10];
    float equinox;
    char raj[16],decj[16];
    float bmaj, bmin, bpa;
    char *start_coord1, *start_coord2;
    char * tmode;
    char *stop_coord1,*stop_coord2;
    double slen;
    char *fdmode;
    float fareq;
    char *calmode;
    double calfreq;
    float caldcyl, calph;
    int calnph;
    long stt_imjd, stt_smjd;
    double stt_offs, stt_lst;
};
struct fromheader {
    double mjd, fract, period;
    long numPulses;
    float freq,dm;
    int numBins,pol,numChan;
    char tid;
};

int open_psrfits(const char * filename);
int put_primary_header();
int put_subint_header();
int put_history_header();
void reportAndExitOnFITSerror(int status);
int readHeader();
int readSubs();
int load_primary_header();
int SuperFileIO();
int FreeThemAll();
int PrintThemAll();


#endif
