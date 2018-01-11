// God, I'm writing C after so many years.
// In the name of stdin, stdout and Dennis Ritchie
// #include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "first.h"
#include "fitsio.h" // This does the magic
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <getopt.h>


// Wow so far so good

/* GLOBAL VARIABLES */
/* Define CFITS file pointer - global- defined in fitsio.h */
fitsfile *fits;
/* records the status all the time */
int sta;
/* File Pointer .. One pointer to rule them all */
FILE * fp;
/* Pointer to pointer to make array of array of header structs */
struct fromheader * fh;
/* Pointer to pointer to make array of arrays of subints */
float * subints;
/* Pointer to PrimaryHDU struct */
struct pHDU * ph;
/* Date time*/
char date_time[24];
/*Here lies the magic constants which at some point will be passed as
  arguments to program.
  TODO: You gotta do it boi!!
*/
int dflag = 0; // Corresponds to Number of Subints `grep '#' | wc -l `
int nflag = 0; // Corresponds to Number of Bins
int bflag = 0;  // Corresponds to Number of Subbands
float bwflag = 0.0; // Corresponds to bandwidth 
char sr[10];

char fn[80],fname[80];
char *file_polyco = "polyco.dat";
/*
    Getting down to business!!!
*/

int load_primary_header(){
    /*
    Some dumb some wise.
    Such are the people and similar are these
    */
    ph = (struct pHDU *)malloc(sizeof(struct pHDU));
    ph->observer = "        ";
    ph->project =  "InPTA";
    ph->telescope = "Ooty";
    ph->ant_x = 1442712.95;
    ph->ant_y = 6087044.73;
    ph->ant_z = 1251052.35;
    ph->fend = "530m Cylindrical Palaboloid";
    ph->ibeam = 1;
    ph->numChan = bflag;
    ph->fd_pol = "LIN";
    ph->fd_had = 1;
    ph->fd_sang = 0.0;
    ph->fd_xy = 0.0;
    ph->bend = "PONDER";
    ph->befcfg = "PONDER_rocks";
    ph->bep = 0;
    ph->bedcc = 0;
    ph->bede = 0.0;
    ph->tcycle = 0.0;
    ph->obsm = "PSR";
    ph->dobs = "2017-05-15T10:06:50.0005";
    ph->obscf = 334.5;
    ph->obsbw = bwflag;
    ph->obsfch = 16; // check
    ph->dm = 100.2;
    ph->pnid = "POINT_ID";
    strcpy(ph->srcn,"YOLO");
    ph->reffrm = "J2000";
    ph->equinox = 2000.0;
    // ph->raj = "00:00:00.0000";
    // ph->decj = "00:00:00.0000";
    ph->bmaj = 0.5;
    ph->bmin = 0.5;
    ph->bpa  = 0.5;
    ph->start_coord1 = "00:00:00.0000";
    ph->start_coord2 = "00:00:00.0000";
    ph->tmode = "TRACK";
    ph->stop_coord1 = "00:00:00.0000";
    ph->stop_coord2 = "00:00:00.0000";
    ph->slen = 12342.123;
    ph->fdmode = "FA";
    ph->fareq = 0.0;
    ph->calmode = "OFF";
    ph->calfreq = 100.0;
    ph->caldcyl = 0.0;
    ph->calph = 0.0;
    ph->calnph = 0;
    ph->stt_imjd = 57872;
    ph->stt_smjd = 23456;
    ph->stt_offs = 0.000123;
    ph->stt_lst = 12.34;
}

int readHeader(){
    /*
    @INPUT struct fromheader * , FILE *
    @OUTPUT 0 if all goes well
    @DOES Puts all header file things in h
    */
    int e4;
    fh = (struct fromheader *)malloc(sizeof(struct fromheader));
    if(fh == NULL){
        return -10;
    }
    e4 = fscanf(fp,"# %lf %lf %lf %ld %f %f %d %c %d \n",
      &(fh->mjd),&(fh->fract),&(fh->period),&(fh->numPulses),
      &(fh->freq),&(fh->dm),&(fh->numBins),&(fh->tid),&(fh->pol));
    fh->numChan = ph->numChan = bflag;
    return e4;
}

int readSubs() {
    /*
    @INPUT FILE Pointer to PROF file && FLOAT Pointer to subs
    @OUTPUT POSITIVE if all goes well
    @DOES Puts the subintegrations in subs
    */
    int e5, iwas, b, c;
    subints = (float*)malloc(fh->numBins*fh->numChan*sizeof(float));
    // HARDCODED FOR SINGLE POLARISATION
    // Simple pointer arithmetic
    for (b=0;b<fh->numBins;b++) {
        e5 = fscanf(fp,"%d",&iwas);
            for (c=0;c<fh->numChan;c++)
                // e5 = fscanf(fp," %f",&subs[b*fh->numChan + c]);
                e5 = fscanf(fp," %f",&subints[c*fh->numBins + b]);
        e5 = fscanf(fp,"\n");
    }
    return e5;
}
int open_psrfits( const char *filename ){
    int bitpix   =  8;
    int naxis   =  0;
    long *naxes;

    sta = 0;
    /* Open file */
    fits_create_file( &fits, filename, &sta );
    /* Write primary header */
    bitpix = 8;
    naxis = 0;
    fits_create_img( fits, bitpix, naxis, naxes, &sta );

    /* Add first keywords */
    fits_update_key( fits, TSTRING, "HDRVER", "5.4", "Header version", &sta);
    fits_update_key(fits, TSTRING, "FITSTYPE","PSRFITS","FITS definition for pulsar data files ", &sta);
    fits_write_date( fits, &sta);
    return( sta );
}

void reportAndExitOnFITSerror(int status) {
  if(status) {
    fits_report_error(stderr, status);
    exit(1);
  }
}

int put_primary_header(){
    fits_update_key(fits, TSTRING, "OBSERVER", ph->observer,"Observer name(s)", &sta);
    fits_update_key(fits, TSTRING, "PROJID", ph->project,"Project name", &sta);
    fits_update_key(fits, TSTRING, "TELESCOP", ph->telescope,"Telescope name", &sta);
    fits_update_key(fits, TDOUBLE, "ANT_X", &(ph->ant_x),"[m] Antenna ITRF X-coordinate (D)", &sta);
    fits_update_key(fits, TDOUBLE, "ANT_Y", &(ph->ant_y),"[m] Antenna ITRF Y-coordinate (D)", &sta);
    fits_update_key(fits, TDOUBLE, "ANT_Z", &(ph->ant_z),"[m] Antenna ITRF Z-coordinate (D)", &sta);
    fits_update_key(fits, TSTRING, "FRONTEND", ph->fend,"Observatory standard name for receiver package", &sta);
    fits_update_key(fits, TINT, "IBEAM",&ph->ibeam,"Beam ID for multibeam systems", &sta);
    fits_update_key(fits, TINT,"NRCVR",&(fh->pol),"Number of receiver polarisation channels", &sta);
    fits_update_key(fits, TSTRING, "FD_POLN", ph->fd_pol,"LIN or CIRC", &sta);
    fits_update_key(fits, TINT,"FD_HAND",&(ph->fd_had),"+/- 1. +1 is LIN:A=X,B=Y, CIRC:A=L,B=R (I) ", &sta);
    fits_update_key(fits, TFLOAT,"FD_SANG",&(ph->fd_sang),"[deg] FA of E vect for equal sig in A&B (E) ", &sta);
    fits_update_key(fits, TFLOAT, "FD_XYPH", &(ph->fd_xy),"[deg] Phase of A* B for injected cal (E) ", &sta);
    fits_update_key(fits, TSTRING, "BACKEND", ph->bend,"Backend ID", &sta);
    fits_update_key(fits, TSTRING, "BECONFIG", ph->befcfg,"Backend configuration file name", &sta);
    fits_update_key(fits, TINT, "BE_PHASE", &(ph->bep), "0/+1/-1 BE cross-phase:0 unknown,+/-1 std/rev ", &sta);
    fits_update_key(fits, TINT, "BE_DCC", &(ph->bedcc), "BE downconversion conjugation corrected", &sta);
    fits_update_key(fits, TFLOAT, "BE_DELAY",&(ph->bede), "[s] Backend propn delay from digitiser input ", &sta);
    fits_update_key(fits, TLONG, "TCYCLE", &(ph->tcycle), "Native cycle time of correlation system",&sta);
    fits_update_key(fits, TSTRING, "OBS_MODE", ph->obsm, "(PSR, CAL, SEARCH) ", &sta);
    fits_update_key(fits, TSTRING, "DATE-OBS", ph->dobs ,"UTC date of observation (YYYY-MM-DDThh:mm:ss) ", &sta);
    double x = (double) ph->obscf;
    fits_update_key(fits, TDOUBLE, "OBSFREQ", &x, "[MHz] Centre frequency for observation",&sta);
    x = (double) ph->obsbw;
    fits_update_key(fits, TDOUBLE, "OBSBW", &x, "[MHz] Bandwidth for observation ",&sta);
    fits_update_key(fits, TINT, "OBSNCHAN", &(ph->obsfch), "Number of frequency channels (original) ",&sta);
    fits_update_key(fits, TFLOAT, "CHAN_DM", &(ph->dm), "[cm-3 pc] DM used for on-line dedispersion ",&sta);
    fits_update_key(fits, TSTRING, "PNT_ID", ph->pnid, "Name or ID for pointing ctr (multibeam feeds) ", &sta);
    fits_update_key(fits, TSTRING, "SRC_NAME", ph->srcn, "Source or scan ID ", &sta);
    fits_update_key(fits, TSTRING, "COORD_MD", ph->reffrm, "Coordinate mode (J2000, GALACTIC, ECLIPTIC) ",&sta);
    fits_update_key(fits, TFLOAT, "EQUINOX", &(ph->equinox), "Equinox of coords (e.g. 2000.0) ",&sta);
    fits_update_key(fits, TSTRING, "RA", ph->raj, "Right ascension (hh:mm:ss.ssss) ",&sta);
    fits_update_key(fits, TSTRING, "DEC", ph->decj, "Declination (-dd:mm:ss.sss) ", &sta);
    fits_update_key(fits, TFLOAT, "BMAJ", &(ph->bmaj), "[deg] Beam major axis length ", &sta);
    fits_update_key(fits, TFLOAT, "BMIN", &(ph->bmin), "[deg] Beam minor axis length ", &sta);
    fits_update_key(fits, TFLOAT, "BPA", &(ph->bpa), "[deg] Beam position angle ", &sta);
    fits_update_key(fits, TSTRING, "STT_CRD1", ph->start_coord1, "Start coord 1 (hh:mm:ss.sss or ddd.ddd) ",&sta);
    fits_update_key(fits, TSTRING, "STT_CRD2", ph->start_coord2, "Start coord 2 (hh:mm:ss.sss or ddd.ddd) ",&sta);
    fits_update_key(fits, TSTRING, "TRK_MODE", ph->tmode, "Track mode (TRACK, SCANGC, SCANLAT) ",&sta);
    fits_update_key(fits, TSTRING, "STP_CRD1", ph->stop_coord1, "Stop coord 1 (hh:mm:ss.sss or ddd.ddd) ",&sta);
    fits_update_key(fits, TSTRING, "STP_CRD2", ph->stop_coord2, "Stop coord 2 (hh:mm:ss.sss or ddd.ddd) ",&sta);
    // fits_update_key(fits, TDOUBLE, "SCANLEN", &(ph->slen), "[s] Requested scan length (E) ",&sta);
    fits_update_key(fits, TSTRING, "FD_MODE", ph->fdmode, "Feed track mode - FA, CPA, SPA, TPA ",&sta);
    fits_update_key(fits, TFLOAT, "FA_REQ", &(ph->fareq), "[deg] Feed/Posn angle requested (E) ",&sta);
    fits_update_key(fits, TSTRING, "CAL_MODE", ph->calmode, "Cal mode (OFF, SYNC, EXT1, EXT2) ",&sta);
    fits_update_key(fits, TDOUBLE, "CAL_FREQ", &(ph->calfreq), "[Hz] Cal modulation frequency (E) ",&sta);
    fits_update_key(fits, TFLOAT, "CAL_DCYC", &(ph->caldcyl), "Cal duty cycle (E) ",&sta);
    fits_update_key(fits, TFLOAT, "CAL_PHS", &(ph->calph), "Cal phase (wrt start time) (E) ",&sta);
    fits_update_key(fits, TINT, "CAL_NPHS", &(ph->calnph), "Number of states in cal pulse (I) ",&sta);
    // printf("In put_primary_header before IO : %ld",ph->stt_imjd);
    fits_update_key(fits, TLONG, "STT_IMJD", &(ph->stt_imjd), "Start MJD (UTC days) (J - long integer)",&sta);
    fits_update_key(fits, TLONG, "STT_SMJD", &(ph->stt_smjd), "[s] Start time (sec past UTC 00h) (J) ",&sta);
    fits_update_key(fits, TDOUBLE, "STT_OFFS", &(ph->stt_offs),"[s] Start time offset (D) ",&sta);
    fits_update_key(fits, TDOUBLE, "STT_LST", &(ph->stt_lst), "[s] Start LST (D) ",&sta);
    fits_update_key(fits, TINT, "NPOL", &(fh->pol), "Number of Polarisations(I)", &sta);
    // fits_get_hdu_num( fits, &last_scanhdr_hdu );
    return(sta);
}

int put_subint_header(){
    // int sta;
    int nscl, ndata;
    nscl = fh->numChan * fh->pol;
    ndata = fh->numBins * nscl;
    int i,j, k;
    float x;
    double dx, dy, dz, ds, dc;

    int ncols=21, col;
    char *ttype[ncols], *tform[ncols], *tunit[ncols];
    long naxes[3];

    char Cstr16[16], Estr16[16], Istr16[16];

    /* Create SUBINT BINTABLE */
    long long nrows = 0; /* naxis2 - Let CFITSIO sort this out */
    /*
    Why keep a standard when you don't follow it?
    Adding a period field in the SUBINT BINTABLE
    */

    ttype[0] = "ISUBINT ";    /* Subint number. If NAXIS=-1, 0 indicates EOD. */
    tform[0] = "1J      ";
    tunit[0] = "";
    ttype[1] = "INDEXVAL";    /* Optionally used if INT_TYPE != TIME */
    tform[1] = "1D      ";
    tunit[1] = "";
    ttype[2] = "TSUBINT ";    /* [s] Length of subintegration */
    tform[2] = "1D      ";
    tunit[2] = "";
    ttype[3] = "OFFS_SUB";    /* [s] Offset from Start UTC of subint centre */
    tform[3] = "1D      ";
    tunit[3] = "";
    ttype[4] = "LST_SUB ";    /* [s] LST at subint centre */
    tform[4] = "1D      ";
    tunit[4] = "";
    ttype[5] = "RA_SUB  ";    /* [turns] RA (J2000) at subint centre */
    tform[5] = "1D      ";
    tunit[5] = "";
    ttype[6] = "DEC_SUB ";    /* [turns] Dec (J2000) at subint centre */
    tform[6] = "1D      ";
    tunit[6] = "";
    ttype[7] = "GLON_SUB";    /* [deg] Gal longitude at subint centre */
    tform[7] = "1D      ";
    tunit[7] = "";
    ttype[8] = "GLAT_SUB";    /* [deg] Gal latitude at subint centre */
    tform[8] = "1D      ";
    tunit[8] = "";
    ttype[9] = "FD_ANG  ";    /* [deg] Feed angle at subint centre */
    tform[9] = "1E      ";
    tunit[9] = "";
    ttype[10] = "POS_ANG ";    /* [deg] Position angle of feed at subint centre */
    tform[10] = "1E      ";
    tunit[10] = "";
    ttype[11] = "PAR_ANG ";    /* [deg] Parallactic angle at subint centre */
    tform[11] = "1E      ";
    tunit[11] = "";
    ttype[12] = "TEL_AZ  ";    /* [deg] Telescope azimuth at subint centre */
    tform[12] = "1E      ";
    tunit[12] = "";
    ttype[13] = "TEL_ZEN ";    /* [deg] Telescope zenith angle at subint centre */
    tform[13] = "1E      ";
    tunit[13] = "";

    sprintf( Cstr16, "%dD", bflag );
    ttype[14] = "DAT_FREQ";
    tform[14] = Cstr16;
    tunit[14] = "";
    ttype[15] = "DAT_WTS ";
    tform[15] = Cstr16;
    tunit[15] = "";

    sprintf( Estr16, "%dE", nscl );
    ttype[16] = "DAT_OFFS";
    tform[16] = Estr16;
    tunit[16] = "";
    ttype[17] = "DAT_SCL ";
    tform[17] = Estr16;
    tunit[17] = "";

    sprintf( Istr16, "%dE", ndata );
    ttype[18] = "DATA    ";
    tform[18] = Istr16;
    tunit[18] = "Jy      ";

    // Y keep a standard if you don't follow it?
    ttype[19] = "PERIOD  ";
    tform[19] = "1D      ";
    tunit[19] = "s       ";
    // Y seriously?
    ttype[20] = "NPOL    ";
    tform[20] = "1I      ";
    tunit[20] = "";

    fits_create_tbl( fits, BINARY_TBL, nrows, ncols, ttype, tform, tunit, "SUBINT", &sta);
    reportAndExitOnFITSerror(sta);

    /* Add dimensions of column 'ncols' = SUBINT Data */
    naxes[0] = fh->numBins;
    naxes[1] = fh->numChan;
    naxes[2] = fh->pol;
    // printf(".....Heere?\n" );
    fits_write_tdim( fits, 19, 3, naxes, &sta );
    reportAndExitOnFITSerror(sta);

    /* Add keywords */
    fits_update_key( fits, TSTRING, "INT_TYPE", "TIME", "Time axis (TIME, BINPHSPERI, BINLNGASC, etc)", &sta);
    fits_update_key( fits, TSTRING, "INT_UNIT", "SEC", "Unit of time axis (SEC, PHS (0-1), DEG)", &sta);
    fits_update_key(fits, TSTRING, "SCALE", "FluxDen","Intensity units (FluxDen/RefFlux/Jansky) ",&sta);
    fits_update_key(fits, TSTRING, "POL_TYPE", "AA", "Polarisation identifier (e.g., AABBCRCI, AA+BB)", &sta);
    fits_update_key(fits, TINT, "NPOL", &(fh->pol), "Nr of polarisations ",&sta);
    x = (float)fh->period / fh->numBins;
    fits_update_key(fits, TFLOAT, "TBIN", &x, "[s] Time per bin or sample ", &sta);
    k = 1;
    fits_update_key(fits, TINT, "NBITS", &k," Number of bits per sample", &sta);
    fits_update_key(fits, TINT, "NBIN", &(fh->numBins), "Nr of bins (PSR/CAL mode; else 1) ", &sta);
    fits_update_key(fits, TINT, "NCHAN", &(fh->numChan), "Number of channels/sub-bands in this file ", &sta);
    dx = (double) ph->obsbw / fh->numChan;
    fits_update_key(fits, TDOUBLE,"CHAN_BW",&dx,"[MHz] Channel/sub-band width ", &sta);
    fits_update_key(fits, TFLOAT, "DM", &(fh->dm), "[cm-3 pc] DM for post-detection dedisperion ", &sta);
    x = 0.0;
    fits_update_key(fits, TFLOAT, "RM", &x, "[rad m-2] RM for post-detection deFaraday", &sta);
    x = 0.0;
    fits_update_key(fits, TFLOAT, "NCHNOFFS", &x, "Channel/sub-band offset for split files ", &sta);
    k = 1; // Search mode
    fits_update_key(fits, TINT, "NSBLK", &k, "Samples/row (SEARCH mode, else 1)", &sta);
    fits_update_key(fits, TINT, "NSTOT", &k, "Total number of samples (SEARCH mode, else 1) ", &sta);
    reportAndExitOnFITSerror(sta);
    // Now focusing on writing the subintegration data.
    int subint_cnt = 1;
    float *binned_weight, *binned_offset, *binned_scale;
    float *binned_freq;

    while(!feof(fp)){
        readSubs();
        nscl = fh->numChan * fh->pol;
        ndata = fh->numBins * nscl;

        binned_weight = (float*)malloc(fh->numChan*sizeof(float));
        binned_offset = (float*)malloc(nscl*sizeof(float));
        binned_scale = (float*)malloc(nscl*sizeof(float));

        for(j = 0; j < fh->numChan;j++){
            binned_weight[j] = 1.0;
        }
        for(j=0;j<nscl;j++){
            binned_offset[j] = 0.0;
            binned_scale[j] = 1.0;
        }
        float dx = ph->obsbw / fh->numChan;
        binned_freq = (float*)malloc(fh->numChan*sizeof(float));
        int ni, nc = fh->numChan/2; // numChan is always a power of 2
        for(ni = 0; ni < fh->numChan; ni++){
            binned_freq[ni] = (float) fh->freq - ((float)ni * dx);
        }
        col =1;
        /* Subint number. If NAXIS=-1, 0 indicates EOD. */
        fits_write_col( fits, TINT, col, subint_cnt, 1, 1, &subint_cnt, &sta );
        col++;
        /* INDEXVAL - Optionally used if INT_TYPE != TIME */
        dx = 0.0;
        fits_write_col( fits, TDOUBLE, col, subint_cnt, 1, 1, &dx, &sta );
        col++;
        /* [s] Length of subint ALAKAZAM */
        double mdx = fh->period * (double)fh->numPulses;
        fits_write_col( fits, TDOUBLE, col, subint_cnt, 1, 1, &mdx, &sta );
        col++;
        /* [s] Offset from Start UTC of subint centre */
        // dx = sum_subint_mid_pt / slen;
        // dx = 86400.0*fh[i]->fract/(2*dx);
        // printf("This is what ? %lf",dx);
		dx = fh->fract/2;
        fits_write_col( fits, TDOUBLE, col, subint_cnt, 1, 1, &(dx), &sta );
        col++;

        /* [s] LST at subint centre */
        // ds = sum_subint_lst_sin / sum_subint_len_secs;
        // dc = sum_subint_lst_cos / sum_subint_len_secs;
        // if( ( dx = ( atan2( ds, dc ) / TwoPi ) ) < 0.0 ) dx += 1.0;
        dx = 86400.0*fh->fract;
        fits_write_col( fits, TDOUBLE, col, subint_cnt, 1, 1, &dx, &sta );
        col++;

        // /* [turns] RA (J2000) at subint centre */
        // ds = sum_subint_ra_sin / sum_subint_len_secs;
        // dc = sum_subint_ra_cos / sum_subint_len_secs;
        // if( ( dx = ( atan2( ds, dc ) / TwoPi ) ) < 0.0 ) dx += 1.0;
        // dx = 1234.1234;
        fits_write_col( fits, TDOUBLE, col, subint_cnt, 1, 1, &dx, &sta );
        col++;

        /* [turns] Dec (J2000) at subint centre */
        // ds = sum_subint_dec_sin / sum_subint_len_secs;
        // dc = sum_subint_dec_cos / sum_subint_len_secs;
        // dx = atan2( ds, dc ) / TwoPi;
        // dx = 1234.1234;
        fits_write_col( fits, TDOUBLE, col, subint_cnt, 1, 1, &dx, &sta );
        col++;

        /* [deg] Gal longitude at subint centre */
        // ds = sum_subint_Glon_sin / sum_subint_len_secs;
        // dc = sum_subint_Glon_cos / sum_subint_len_secs;
        // if( ( dx = ( atan2( ds, dc ) / TwoPi ) ) < 0.0 ) dx += 1.0;
        // dx = 360.0;
        fits_write_col( fits, TDOUBLE, col, subint_cnt, 1, 1, &dx, &sta );
        col++;

        /* [deg] Gal latitude at subint centre */
        // ds = sum_subint_Glat_sin / sum_subint_len_secs;
        // dc = sum_subint_Glat_cos / sum_subint_len_secs;
        // dx = atan2( ds, dc ) * 360.0 / TwoPi;
        // dx = 360.0;
        fits_write_col( fits, TDOUBLE, col, subint_cnt, 1, 1, &dx, &sta );
        col++;

        /* [deg] Feed angle at subint centre */
        // ds = sum_subint_fa_sin / sum_subint_len_secs;
        // dc = sum_subint_fa_cos / sum_subint_len_secs;
        // dx = atan2( ds, dc ) * 360.0 / TwoPi;
        dx = 360.0;
        x = (float) dx;
        fits_write_col( fits, TFLOAT, col, subint_cnt, 1, 1, &x, &sta );
        col++;

        /* [deg] Parallactic angle at subint centre */
        // ds = sum_subint_pa_sin / sum_subint_len_secs;
        // dc = sum_subint_pa_cos / sum_subint_len_secs;
        // dy = atan2( ds, dc ) * 360.0 / TwoPi;

        /* [deg] Position angle of feed at subint centre */
        // dx = dx + dy;
        // if( dx > 180.0 ) dx -= 360.0;
        // if( dx < 180.0 ) dx += 360.0;
        // x = (float) dx;
        fits_write_col( fits, TFLOAT, col, subint_cnt, 1, 1, &x, &sta );
        col++;

        /* [deg] Parallactic angle at subint centre */
        // x = (float) d;
        fits_write_col( fits, TFLOAT, col, subint_cnt, 1, 1, &x, &sta );
        col++;

        /* [deg] Telescope azimuth at subint centre */
        // ds = sum_subint_az_sin / sum_subint_len_secs;
        // dc = sum_subint_az_cos / sum_subint_len_secs;
        // if( ( dx = ( atan2( ds, dc ) / TwoPi ) ) < 0.0 ) dx += 1.0;
        // dx = 360.0;
        x = (float) dx;
        fits_write_col( fits, TFLOAT, col, subint_cnt, 1, 1, &x, &sta );
        col++;

        /* [deg] Telescope zenith angle at subint centre */
        // ds = sum_subint_el_sin / sum_subint_len_secs;
        // dc = sum_subint_el_cos / sum_subint_len_secs;
        // dx = 90.0 - ( atan2( ds, dc ) * 360.0 / TwoPi );
        // x = (float) dx;
        // x  = 360.0;
        fits_write_col( fits, TFLOAT, col, subint_cnt, 1, 1, &x, &sta );
        col++;
        /* Centre freq. for each channel - NCHAN floats */
        int e8;
        fits_write_col( fits, TFLOAT, col, subint_cnt, 1, fh->numChan, binned_freq, &sta );
        col++;
        /* Weights for each channel -  NCHAN floats */
        fits_write_col( fits, TFLOAT, col, subint_cnt, 1, fh->numChan, binned_weight, &sta );
        col++;
        /* Data offset for each channel - NCHAN*NPOL floats */
        fits_write_col( fits, TFLOAT, col, subint_cnt, 1, nscl, binned_offset, &sta );
        col++;
        /* Data scale factor for each channel - NCHAN*NPOL floats */
        fits_write_col( fits, TFLOAT, col, subint_cnt, 1, nscl, binned_scale, &sta );
        col++;
        /* Subint data table - Dimensions of data table = (NBIN,NCHAN,NPOL) */
        fits_write_col( fits, TFLOAT, col, subint_cnt, 1, ndata, subints, &sta );
        col++;
        reportAndExitOnFITSerror(sta);
        /* Why keep a standard when you have to break it? */
        fits_write_col(fits, TDOUBLE, col, subint_cnt, 1, 1, &(fh->period), &sta);
        col++;
        /* Why keep a standard when you have to break it? */
        fits_write_col(fits, TINT, col, subint_cnt, 1, 1, &(fh->pol), &sta);
        col++;
        subint_cnt++;
        free(binned_weight);free(binned_offset);free(binned_scale);
        free(binned_freq);
        free(fh);
        free(subints);
        readHeader(); // ensure fh->numChan is loaded
        dflag++;
        /* Finished SUBINT */
    }
    reportAndExitOnFITSerror(sta);
    return(sta);
}

int put_history_header(){
    /*
    Missed the following here :
    fd_corr, be_corr, rm_model, aux_rm_c, dm_model, aux_dm_c
    When writing TSTRING to fits_write_col you gotta send pointer to char*.
    */
    int col = 1;
    int n = 0;
    char *PHtype[20] = {
      "DATE_PRO","PROC_CMD","POL_TYPE","NPOL    ","NBIN    ","NBIN_PRD","TBIN    ","CTR_FREQ",
      "NCHAN   ","CHAN_BW ","PAR_CORR","RM_CORR ","DEDISP  ","DDS_MTHD","SC_MTHD ","CAL_MTHD",
      "CAL_FILE","RFI_MTHD" // ,"SCALE","NSUB"
    };
    char *PHform[20] = {
      "24A     ","80A     ","8A      ","1I      ","1I      ","1I      ","1D      ","1D      ",
      "1I      ","1D      ","1I      ","1I      ","1I      ","32A     ","32A     ","32A     ",
      "32A     ","32A     " // ,"8A      ","1I      "
    };
    char *PHunit[20] = {
      "        ","        ","        ","        ","        ","        ","s       ","MHz     ",
      "        ","MHz     ","        ","        ","        ","        ","        ","        ",
      "        ","        " // ,"        ","        "
    };
    int nrows = 0;  /* naxis2 - Let CFITSIO sort this out */
    int ncols = 18; /* tfields */
    fits_create_tbl( fits, BINARY_TBL, nrows, ncols, PHtype, PHform, PHunit, "HISTORY", &sta);
    /* Processing date and time (YYYY-MM-DDThh:mm:ss UTC) */
    char * mystring;
    mystring = date_time;
    fits_write_col( fits, TSTRING, 1, 1, 1, 1, &mystring, &sta );
    /* Processing program and command */
    mystring = "WBCORR";
    fits_write_col( fits, TSTRING, 2, 1, 1, 1, &mystring, &sta );
    /* Polarisation identifier */
    mystring = "AA";
    fits_write_col( fits, TSTRING, 3, 1, 1, 1, &mystring, &sta );
    /* Nr of pols*/
    fits_write_col( fits, TSHORT, 4, 1, 1, 1, &(fh->pol), &sta );
    /* Nr of bins per product (0 for SEARCH mode) */
    fits_write_col( fits, TSHORT, 5, 1, 1, 1, &(fh->numBins), &sta );
    /* Nr of bins per period */
    // TODO: Why is same thing written again? check
    fits_write_col( fits, TSHORT, 6, 1, 1, 1, &(fh->numBins), &sta );
    /* Bin time */
    double dx = (double)fh->period / fh->numBins;
    fits_write_col( fits, TDOUBLE, 7, 1, 1, 1, &dx, &sta );
    /* Centre freq. */
    // dx = (double)(fh->freq);
    dx = (double) ph->obscf;
    fits_write_col( fits, TDOUBLE, 8, 1, 1, 1, &dx, &sta );
    /* Number of channels */
    fits_write_col( fits, TSHORT, 9, 1, 1, 1, &(fh->numChan), &sta );
    /* Channel bandwidth */
    dx = - (double)ph->obsbw / fh->numChan;
    // printf("Channel Bandwidth %lf\n", dx);
    fits_write_col( fits, TDOUBLE, 10, 1, 1, 1, &dx, &sta );
    int e3 = 0;
    /* Parallactic angle correction applied */
    fits_write_col( fits, TSHORT, 11, 1, 1, 1, &e3, &sta );
    /* RM correction applied */
    fits_write_col( fits, TSHORT, 12, 1, 1, 1, &e3, &sta );
    /* Data dedispersed */
    e3 = 1;
    fits_write_col( fits, TSHORT, 13, 1, 1, 1, &e3, &sta );
    /* Dedispersion method */
    mystring = "SIGPROC";
    fits_write_col( fits, TSTRING, 14, 1, 1, 1, &mystring, &sta );
    /* Scattered power correction method */
    mystring = "SCATTER";
    fits_write_col( fits, TSTRING, 15, 1, 1, 1, &mystring, &sta );
    /* Calibration method */
    mystring = "NONE";
    fits_write_col( fits, TSTRING, 16, 1, 1, 1, &mystring, &sta );
    /* Name of calibration file */
    mystring = "NONE";
    fits_write_col( fits, TSTRING, 17, 1, 1, 1, &mystring, &sta );
    /* RFI excision method */
    mystring = "ZAP";
    fits_write_col( fits, TSTRING, 18, 1, 1, 1, &mystring, &sta );
    // /* SCALE */
    // mystring = "FluxDen";
    // fits_write_col(fits, TSTRING, 19, 1, 1, 1, &mystring, &sta);
    // /* NSUB */
    // printf("NSUBS zero ? %d \n", fh->numBins);
    // fits_write_col(fits, TINT, 20, 1,1,1, &(fh->numBins), &sta);
    return(sta);
}

int
main(int argc, char *argv[]) {
    FILE * fpar;
    char fparn[16];
    char praj[16],pdecj[16];
    int fio = 0,opt,n;
    char *ofn;
    fn[0] = '\0';
    fname[0] = '\0';
    if(argc < 2){
        printf("\t -------------PSRFITS Generator-------------\n");
        printf("Following are the options you'll have to pass.\n");
        printf("\t\t-d  --> Number of Subintegrations that will be there.\n");
        printf("\t\t\t\tIf not given, will be determined on the fly.\n" );
        printf("\t\t-b  --> Number of Sub-bands.(Same as what was given to dedisperse)\n");
        printf("\t\t-n  --> Number of bins.(Same as what was given to foldT2)\n");
        printf("\t\t-o  --> Name of the output FITS file.(Should be unique)\n");
        printf("\t\t-i  --> Name of the input file.(default is stdin)\n");
        printf("\t\t-s  --> Source Name\n");
        printf("\t\t-u  --> Right Ascesion of source\n");
        printf("\t\t-v  --> Declination of source\n");
        printf("\t\t-f  --> Config file\n");
        printf("\t\t-w  --> Bandwidth\n");
        printf("\t\t ---------------------\n");
        printf("Written by Suryrao Bethapudi(ep14btech11008@iith.ac.in)\n");
        printf("\t\t ---------------------\n");
        exit(0);
    }
    while( (opt = getopt(argc,argv,"d:b:n:o:i:f:u:v:s:w:")) != -1 ){
        switch(opt) {
		    case 'w' :
		    	bwflag = atof(optarg);
		    	if(bwflag <= 1){
		    	    fprintf(stderr,"Argument passed to -w is erroneous : %f.\n",bwflag);
		    	    exit(-3);
		    	}
            case 'd' :
                dflag = atoi(optarg);
                if(dflag <=0){
                    fprintf(stderr,"Argument passed to -d is erroneous : %d.\n",dflag);
                    exit(-3);
                }
                break;
            case 'b' :
                bflag = atoi(optarg);
                if(bflag <=0){
                    fprintf(stderr,"Argument passed to -b is erroneous : %d.\n",bflag);
                    exit(-3);
                }
                break;
            case 'n' :
                nflag = atoi(optarg);
                if(nflag <=0){
                    fprintf(stderr,"Argument passed to -n is erroneous : %d.\n",nflag);
                    exit(-3);
                }
                break;
            case 'o' :
                // fname = optarg;
                strcpy(fname,optarg);
                ofn = strcat(fname,".fits");
                if(access(fname, F_OK) != -1){
                    fprintf(stderr,"FITS file with that name already exists: %s\n",ofn);
                    fprintf(stderr,"Consider deleting it.\n");
                    exit(-3);
                }
                break;
            case 'i' :
                // fn = optarg;
                strcpy(fn,optarg);
                if(access(fn, F_OK) == -1){
                    fprintf(stderr,"No file with that name exists: %s\n",fn);
                    exit(-3);
                }
                fio = 1;break;
            case 's' :
                strcpy(sr,optarg);
                break;
            case 'f':
                strcpy(fparn,optarg);
                if(access(fparn,F_OK) == -1){
                    fprintf(stderr,"No file with that name exists: %s\n",fparn);
                    exit(-3);
                }
                fpar = fopen(fparn,"r");
                // FPAR IO
                fscanf(fpar,"%s\n",sr); // source name
                fscanf(fpar,"%s\n",praj); // RAJ
                fscanf(fpar,"%s\n",pdecj); // DECJ
                fscanf(fpar,"%d\n",&bflag); // subbands
                fscanf(fpar,"%d\n",&nflag); // numbins
                fclose(fpar);
                break;
            case 'u':
                strcpy(praj,optarg);
                break;
            case 'v':
                strcpy(pdecj,optarg);
                break;
            default:
                fprintf(stderr,"Option not recognized.");
                exit(-2);
        }
    }
    load_primary_header();
    // if(fio == 0){
    //     // STDIN input
    //     SuperPipeIO();
    // }
    // else{
    //     // File Input
    //     fp = fopen(fn,"r");
    //     SuperFileIO();
    // }
    fp = fopen(fn,"r");
    readHeader(); // fp
    // fh = fh[0];
    // ph->slen = (double)dflag * fh->period * (double)fh->numPulses;
    ph->stt_imjd = (long)fh->mjd;
    ph->stt_smjd = (long)fh->fract;
    ph->stt_offs = (fh->fract - (double)ph->stt_smjd);
    ph->dm = fh->dm;
    fits_get_system_time( date_time, &n, &sta );
    ph->dobs = date_time;
    ph->obsfch = fh->numChan = bflag;
    strcpy(ph->srcn,sr);
    strcpy(ph->raj,praj);
    strcpy(ph->decj,pdecj);
    open_psrfits(ofn);
    reportAndExitOnFITSerror(sta);
    put_primary_header();
    reportAndExitOnFITSerror(sta);
    put_history_header();
    reportAndExitOnFITSerror(sta);
    put_subint_header();
    reportAndExitOnFITSerror(sta);
    // To write that scanlen.
    int last_scanhdr_hdu;
    fits_get_hdu_num( fits, &last_scanhdr_hdu ); // Getting last scanned HDU
    fits_movabs_hdu( fits, 1, NULL, &sta ); // moving to primayHDU
    ph->slen = (double)dflag * fh->period * (double)fh->numPulses;
    fits_update_key(fits, TDOUBLE, "SCANLEN", &(ph->slen), "[s] Requested scan length (E) ",&sta);
    reportAndExitOnFITSerror(sta);
    fits_movabs_hdu( fits, last_scanhdr_hdu, NULL, &sta ); // Re-moving back to last scanned HDU
    fits_close_file(fits,&sta);
return sta;
}
