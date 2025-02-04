MODULE rings_simple
    USE precision
    USE data_input
    USE neighbor_finder, only: neigh_list
    USE stdlib_array
    IMPLICIT NONE
    TYPE ring
        integer :: l
        integer :: element(20)
        integer :: sorted(20)
    END TYPE

    integer :: ring_cap, ring_size  ! Main ring list capacity and current size.
    real(dp) :: tstart, tcheck, tneilist, tfindring, tcheckrepi, tcheckpr, taddring
    integer :: crude_ring_num

    type(ring) :: emptyring

CONTAINS

void dist(real x1, real y1, real z1, real x2, real y2, real z2, real *box, real *pbc, real *r) {
  real xij,yij,zij;
  xij=(x1-x2);
  if (pbc[0]==1.0) { if (xij+xij>box[0]) xij-=box[0]; else if (xij+xij<-box[0]) xij+=box[0]; }
  yij=(y1-y2);
  if (pbc[1]==1.0) { if (yij+yij>box[1]) yij-=box[1]; else if (yij+yij<-box[1]) yij+=box[1]; }
  zij=(z1-z2);
  if (pbc[2]==1.0) { if (zij+zij>box[2]) zij-=box[2]; else if (zij+zij<-box[2]) zij+=box[2]; }
  *r=sqrt(xij*xij+yij*yij+zij*zij);
}

subroutine dijk
  implicit none
  do t = 1, nat
    lvldist(t) = lvlreq+1
  end do
  lvldist(nodsrc) = 0
  queue(0) = nodsrc
  dat = 0

  do while (quebgn<=quend)
    nodcrt = queue(quebgn)
    lvlprb = lvldist(nodcrt)+1
    do lnkscrt = 0, lnks(nodcrt)
      nodprb = nodlnkd(nodcrt*mxlinks+lnkscrt)
      if (lvldist(nodprb) > lvlprb) then
        lvldist(nodprb) = lvlprb
        if (lvlprb < lvlreq) then
          queue(quend)=nodprb
          quend = quend+1
          dat = dat+1
        end if
      end if
    end do
  end do

end subroutine

! void dijk(int nodsrc,int lvlreq, int *lnks,int *nodlnkd,int *lvldist,int nat) {
!   int t;
!   int quebgn,quend,nodcrt;
!   int lnkscrt,nodprb,lvlprb;
!
!   for (t=0;t<nat;t++) {lvldist[t]=lvlreq+2;}
!   lvldist[nodsrc]=0;
!   queue[0]=nodsrc;
!   dat=0;
!
!   for (quebgn=0,quend=1;quebgn<=quend;quebgn++) {
!     nodcrt=queue[quebgn];
!     lvlprb=lvldist[nodcrt]+1;
!     for (lnkscrt=0;lnkscrt<lnks[nodcrt];lnkscrt++) {
!       nodprb=nodlnkd[nodcrt*mxlinks+lnkscrt];
!       if (lvldist[nodprb]>lvlprb) {
!         lvldist[nodprb]=lvlprb;
!         if (lvlprb<lvlreq) {
!           queue[quend]=nodprb; // MB: increment 'quend' AFTER adding node to
!           quend++;             // queue, otherwise 2nd node will always be
!           dat++;               // node 0 (cf. Yuan 2002, arrays start from 1)
!         }
!       }
!     }
!   }
!
! }

subroutine srt_rec
  implicit none
  do lnkcrt = 1, lnks(nodcrt)
    nodprb = nodlnkd(nodcrt*mxlinks+lnkcrt)
    lvlprb = lvldist(nodprb)
    if (lvlprb == 0) then
      do lvl = 1, lvlprim+1
        srtpth
      end do
    end if
  end do
end subroutine

int srtpthx[2048];
void srt_rec(int nodcrt,int lvlcrt,int lvlprim){
  //  int srtpthx[mxlevel/2];
  int lnkcrt,nodprb,lvlprb,lvl;
  int t;

  for (lnkcrt=0;lnkcrt<lnks[nodcrt];lnkcrt++) {
    nodprb=nodlnkd[nodcrt*mxlinks+lnkcrt];
    lvlprb=lvldist[nodprb];
    if (lvlprb==0) {
      for (lvl=1;lvl<=lvlprim+1;lvl++) {
	srtpth[((pths)*(mxlevel/2))+lvl]=srtpthx[lvl];
      }
      pths++;
    } else if (lvlprb==(lvlcrt-1)) {
      srtpthx[lvlprb]=nodprb;
      srt_rec(nodprb,lvlprb,lvlprim);
    }
  }
}


void prime_ring(int iodd,int lvlprim,int rngs,int rat) {

  int i,j,k,u;
  struct lit limit;
  int goal_found,irng,pth1,pth2,nodchk,lvlmax,lvlchk;
  int nodmid,irgx,p1x,p2x,t,rlev;

  goal_found=0;

  for (irng=0;irng<rngs;irng++) {
    pth1=querng[irng];
    pth2=querng[irng+mxpaths*mxpaths/2];
    if (pth1>=0) {
      for (lvlmax=lvlprim;lvlmax<=(lvlprim+iodd);lvlmax++) {
	for (lvlchk=1;lvlchk<=(lvlmax-1);lvlchk++) {           // onko herra bugi täällä??
	//for (lvlchk=0;lvlchk<(lvlmax-1);lvlchk++) {           // Which one should it be????
	  nodchk=srtpth[pth1*(mxlevel/2)+lvlchk];
          nodmid=srtpth[pth2*(mxlevel/2)+lvlmax-lvlchk];

	  for (i=1;i<=4;i++) limit.i[i][1]=lvlref[nodmid*4+i]+lvlprim-1;
	  for (i=1;i<=4;i++) limit.i[i][2]=lvlref[nodmid*4+i]-lvlprim+1;
          limit.i[5][1]=lvldist[nodmid]+lvlprim-1;
          limit.i[5][2]=lvldist[nodmid]-lvlprim+1;
          pair_search(&nodchk,&nodmid,&limit,&goal_found);

          if (goal_found==1) {
            goal_found=0;
            for (irgx=irng+1;irgx<rngs;irgx++) {
	      p1x=querng[irgx];
	      p2x=querng[irgx+mxpaths*mxpaths/2];
	      if (p1x>=0) {
		if ((srtpth[p1x*(mxlevel/2)+lvlchk]==nodchk) && srtpth[p2x*(mxlevel/2)+lvlmax-lvlchk]==nodmid) querng[irgx]=-1;
		if ((srtpth[p2x*(mxlevel/2)+lvlchk]==nodchk) && srtpth[p1x*(mxlevel/2)+lvlmax-lvlchk]==nodmid) querng[irgx]=-1;
	      }
            }
	    goto nring;
          } // goal found
	}
      } // lvlmax loop
      rlev=2*lvlprim+iodd;
      ringstat[rlev]=ringstat[rlev]+1;
      j=srtpth[pth1*(mxlevel/2)+1];
      k=srtpth[pth2*(mxlevel/2)+1];
      for (u=0;u<lnks[rat];u++) {
	if (nodlnkd[rat*mxlinks+u]==j)
	  if (sring[rat*mxlinks+u]>rlev) sring[rat*mxlinks+u]=rlev;
        if (nodlnkd[rat*mxlinks+u]==k)
	  if (sring[rat*mxlinks+u]>rlev) sring[rat*mxlinks+u]=rlev;
      }
    } // pth1 > 0
  nring:;
  }
}


void pair_search(int *nodcrt,int *nodgoal,struct lit *limit, int *goal_found) {
  int lnkscrt,nodprb,iref;
  struct lit lmtx;

  if (*nodcrt==*nodgoal) {
    *goal_found=1;
  } else {
      for (lnkscrt=0;lnkscrt<lnks[*nodcrt];lnkscrt++) {
	nodprb=nodlnkd[*nodcrt*mxlinks+lnkscrt];
	for (iref=1;iref<=4;iref++) {
	  if (lvlref[nodprb*4+iref]>=limit->i[iref][1]||
	     lvlref[nodprb*4+iref]<=limit->i[iref][2]) goto nsear;
       	} // iref loop
        if (lvldist[nodprb]>=limit->i[5][1] || lvldist[nodprb]<=limit->i[5][2]) goto nsear;
        for (iref=1;iref<=5;iref++) {
	  lmtx.i[iref][1]=limit->i[iref][1]-1;
          lmtx.i[iref][2]=limit->i[iref][2]+1;
	}
	pair_search(&nodprb,nodgoal,&lmtx,goal_found);
        if (*goal_found==1) return;
      nsear:;
      } // lnkscrt loop
    } // if nodcrt
}



void ringanal(
	      real *x0,
	      int *itype,
	      int *ident,
	      real *box,
	      real *pbc,
	      real cutR[MAXTYPE][MAXTYPE],
	      real cutS[MAXTYPE][MAXTYPE],
	      int nat,
	      real a,
	      real dr,
	      real time,
	      real printshort,
	      logical printangles,
	      logical writenbonds,
	      int printnot,
	      real printeffectivea,
	      real effectivearmax,
	      real effectiveafact
	      )

/*
  Remember that due to the nature of arrays in C and Fortran,
  any index in C is one less than that in Fortran !
*/

{


  static int analtime=0;   /* Counter of how many time we've been here */

  /* Array which returns neighbour info */
  struct ngbr ngbrs[MAXNGBRS];
  int nngbrstat[MAXNGBRS];
  int typenngbrstat[MAXTYPE][MAXNGBRS];
  int typenngbrstattypestat[MAXTYPE][MAXNGBRS][MAXTYPE];
  int numberoftypetypepairs[MAXTYPE][MAXNGBRS][MAXNGBRS];

  int oneatomngbrtypestat[MAXTYPE];

  char file[MAXBUF],bdfile[MAXBUF];

  int jj,kk,nngot,nrave,it,jt,u,v;
  real r,boxmin,rave,r110ave,r1,r2,r3,theta;

  /* Neighbour parameters */
  static real rnmax;
  int npairs;         /* Total number of neighbours */

  int drstat[1000000],ndr;

  /* Cutoff function params */
  real nngbr,nngbrmax,fc;

  int i,i3,j,k;
  int t1,t2,t3,t4,t5;

  static logical firsttime=True;

  // ***************************************************************
  // RING
  // ***************************************************************

  int nodcrt,lvlcrt,lvlprim;

  rnmax=cutS[0][0];
  for (i=0;i<MAXTYPE;i++) for (j=0;j<MAXTYPE;j++) {
      if (cutS[i][j] > rnmax) rnmax=cutS[i][j];
    }
  printf("Starting ring analysis at time %g for %d atoms. rcut %g\n",time,nat,rnmax);

  /* Initialize neighbour list calc. */
  i=-1;
  getngbrs(i,x0,box,pbc,nat,rnmax,ngbrs,&nngot,ident);
  npairs=0;

  /* Initialize statistics */
  for (i=0;i<MAXNGBRS;i++) {
    nngbrstat[i]=0;
  }
  for (i=0;i<MAXTYPE;i++) {
    oneatomngbrtypestat[i]=0;
    for (j=0;j<MAXNGBRS;j++) {
      typenngbrstat[i][j]=0;
      for (k=0;k<MAXTYPE;k++) {
	typenngbrstattypestat[i][j][k]=0;
      }
      for (k=0;k<MAXNGBRS;k++) {
	numberoftypetypepairs[i][j][k]=0;
      }
    }
  }
  boxmin=box[0];
  if (box[1]<boxmin) boxmin=box[1];
  if (box[2]<boxmin) boxmin=box[2];
  ndr=boxmin/dr+1; if (ndr>1000000) { printf("Increase drstat size %d",ndr); exit(0); }
  for (i=0;i<ndr;i++) drstat[i]=0;
  nngbrmax=0.0;


  // Create neighbour list

  printf("Creating neighbour list\n");
  for (i=0;i<nat;i++) {
    getngbrs(i,x0,box,pbc,nat,rnmax,ngbrs,&nngot,ident);
    lnks[i]=nngot;
    t1=0;
    t3=itype[i];
    for (jj=0;jj<nngot;jj++)  {
      t2=itype[ngbrs[jj].i];
      if (cutS[t3][t2]>=ngbrs[jj].r) {
	nodlnkd[i*mxlinks+t1]=ngbrs[jj].i;
	t1++;
      }
    }
    lnks[i]=t1;
    npairs+=nngot;
  }


  for (i=3;i<mxlevel;i++) ringstat[i]=0;
  for (i=0;i<(nat*mxlinks);i++) sring[i]=mxlevel+1;

  printf("Creating 4-point guides\n");
  dijk(nat/4*0,nat,lnks,nodlnkd,lvldist,nat);
  for (i=0;i<nat;i++) lvlref[i*4+1]=lvldist[i];
  dijk(nat/4*1,nat,lnks,nodlnkd,lvldist,nat);
  for (i=0;i<nat;i++) lvlref[i*4+2]=lvldist[i];
  dijk(nat/4*2,nat,lnks,nodlnkd,lvldist,nat);
  for (i=0;i<nat;i++) lvlref[i*4+3]=lvldist[i];
  dijk(nat/4*3,nat,lnks,nodlnkd,lvldist,nat);
  for (i=0;i<nat;i++) lvlref[i*4+4]=lvldist[i];

  pat=0;

  // ring statistics

  for (i=0;i<nat;i++){
    if ((i%1000)==0) printf("%d",i/1000);
    if ((i%200)==0) {printf(".");fflush(stdout);}
    // ==AK== : Statistics only for rings starting at the desired atom type
    if (!start_any && strcmp(start_name,name[itype[i]])!=0) continue;

    // create shortest distance map
    dijk(i,(mxlevel/2)+1,lnks,nodlnkd,lvldist,nat);

    for (t5=0;t5<=dat;t5++) {
      jj=queue[t5];
      if (lvldist[jj]<(mxlevel/2)) {
	v=0;
	for (k=0;k<lnks[jj];k++)
	  if ((lvldist[jj]-1)==lvldist[nodlnkd[jj*mxlinks+k]]) v++;
	if (v>1) {
	  nodcrt=jj;
	  lvlcrt=lvldist[jj];
	  lvlprim=lvlcrt;
	  pths=0;
	  srtpthx[lvlprim]=nodcrt;
	  srt_rec(nodcrt,lvlcrt,lvlprim);
	  if (pat<pths) pat=pths;

	  t4=0;
	  for (t1=0;t1<pths;t1++)
	    for (t2=t1+1;t2<pths;t2++) {
	      t3=0;
	      for (u=1;u<lvlprim;u++)
		if (srtpth[((t1)*mxlevel/2)+u]==srtpth[((t2)*mxlevel/2)+u]) {t3=1;u=lvlprim+1;}
	      if (t3==0) {
		querng[t4]=t1;
		querng[t4+mxpaths*mxpaths/2]=t2;
		t4++;
		if (t4>(mxpaths*mxpaths/2)) {printf("mxpaths overflow:  %d >  %d\n",t4,mxpaths*mxpaths/2);exit(-1);}
	      }
	    }
	  if (t4>0) {
	    prime_ring(0,lvlprim,t4,i);
	  }
	}

	for (k=0;k<lnks[jj];k++)
	  if (lvldist[jj]==lvldist[u=nodlnkd[jj*mxlinks+k]])
	    if (u<jj) {
	      nodcrt=jj;
	      lvlcrt=lvldist[jj];
	      lvlprim=lvlcrt;
	      pths=0;
	      srtpthx[lvlprim]=nodcrt;
	      srt_rec(nodcrt,lvlcrt,lvlprim);
	      j=pths;
	      nodcrt=u;
	      lvlcrt=lvldist[jj];
	      lvlprim=lvlcrt;
	      srtpthx[lvlprim]=nodcrt;
	      srt_rec(nodcrt,lvlcrt,lvlprim);
	      if (pat<pths) pat=pths;
	      t4=0;
	      for (t1=0;t1<j;t1++)
		for (t2=j;t2<pths;t2++) {
		  t3=0;
		  for (u=1;u<lvlprim;u++)
		    if ( srtpth[((t1)*mxlevel/2)+u]==srtpth[((t2)*mxlevel/2)+u]) {t3=1;u=lvlprim+1;}
		  if (t3==0) {
		    querng[t4]=t1;
		    querng[t4+mxpaths*mxpaths/2]=t2;
		    t4++;
		    if (t4>(mxpaths*mxpaths/2)) {printf("mxpaths overflow\n");exit(-1);}
		  }
		} // t2,t1
	      if (t4>0) {
		prime_ring(1,lvlprim,t4,i);
	      }
	    }
      }
    } // t5: 0->dat

    //printf("querng %d : ",i);
    //for (int ii=0;ii<(mxpaths+1)*(mxpaths+1);ii++) printf("%d ",querng[ii]);
    //printf("\n");

  } // i: 0->nat

  printf("\n");
  analtime++;
  printf("Ring statistics : \n");
  for (i=3;i<mxlevel;i++)
    printf("Size %d rings %d\n",i,ringstat[i]/i);

  for (i=0;i<mxlevel;i++)
    ringstat[i]=0;
  for (i=0;i<nat;i++)
    for (t1=0;t1<lnks[i];t1++)
      ringstat[sring[i*mxlinks+t1]]++;

  printf("\nSmallest rings for bonds : \n");
  for (i=3;i<mxlevel;i++)
    printf("SSize %d bonds %d\n",i,ringstat[i]/2);

  fflush(stdout); fflush(stderr);
  firsttime=False;

}

END MODULE rings_simple
