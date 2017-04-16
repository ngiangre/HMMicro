#include "common.h"
#include "geneGraph.h"
#include "genePred.h"
#include "altGraphX.h"
#include "ggMrnaAli.h"
#include "psl.h"
#include "dnaseq.h"
#include "hdb.h"
#include "jksql.h"
#include "cheapcgi.h"
#include "bed.h"



char *agTemplate[] = {"chr22","15615852","15616760","NA","0","-","4","0,3,1,2,","15615854,15616296,15616606,15616755,","3","0,1,2,","1,2,3,","{59,{58,57,47,48,49,50,51,52,53,54,55,56,46,45,44,43,41,42,40,39,38,37,31,32,33,34,35,36,29,30,27,28,24,25,26,23,22,21,20,19,18,17,16,15,14,7,8,9,10,11,12,0,1,2,3,4,5,6,13,},},{59,{58,57,56,55,54,53,52,51,50,49,48,47,46,45,44,43,42,41,40,39,38,37,36,35,34,33,32,31,30,29,28,27,26,25,24,23,22,21,20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0,},},{59,{56,47,55,45,43,39,41,42,38,58,37,36,35,34,48,30,26,31,32,33,40,44,49,51,52,53,57,25,24,23,28,22,20,19,18,27,29,46,17,16,15,14,13,12,11,10,8,50,7,21,6,5,9,3,4,2,1,0,54,},},","1,0,0,","59","BE042632,AW264081,AI830079,AI800406,AI092874,AI016636,AA780238,AW512212,AI971040,AI422048,AI375628,AA778882,H98673,AW674617,AI802847,W78104,F21278,AA485621,C21211,AA115160,AA649623,AI282867,AA701239,F08896,H09581,H09130,BG611378,AW023908,AA977905,AW513448,AI654250,AW027289,AI143197,AI050838,Z38903,F02225,AA022508,AA583258,W77841,AA023025,AI707982,AA715555,AA714920,H92141,AA446668,AA256418,D25564,BE147539,AW813964,AW390758,AA730569,BF839675,BF839672,BF839629,BF839614,BF350389,BE706229,BF839621,AW370225,","3226599,4019974,3616477,3218583,0,0,0,3217336,0,616421,616421,0,84355,0,3637809,0,2278530,212912,0,0,383501,3226599,0,2258019,0,0,127778,231423,397250,4689264,3224495,0,0,0,2258019,2258019,0,97962,0,0,0,97646,97646,1410,0,391011,0,0,0,0,97646,0,0,0,0,0,0,0,0,","97363,4019972,3616475,97309,89851,97123,97061,97313,202617,97075,97075,84361,84353,26677,3637807,89851,2284374,135559,2196729,3454235,97142,97363,84361,2258017,82142,82142,59401,1817081,397248,4689262,97570,3706790,89851,97169,2258017,2258017,89851,6245700,89851,89851,3522683,135621,135621,23404,97061,97245,2223959,523195,586069,586069,135621,523195,523195,523195,523195,4922449,523195,523195,4911150,"};


void usage() 
{
errAbort("testAltGraphXMem - quick little program to test the memory allocation and\n"
	 "freeing in altGraphX library. Must be stopped manually.\n"
	 "usage:\n"
	 "\ttestAltGraphXMem test");
}

int main(int argc, char *argv[])
{
struct altGraphX *ag = NULL;
if(argc != 2)
    usage();
while(1)
    {
    ag = altGraphXLoad(agTemplate);
    altGraphXFree(&ag);
    }
}
