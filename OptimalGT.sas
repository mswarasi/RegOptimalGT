
/================================================*/
/* Appendix J: SAS PROC IML code for analyzing   */
/* the IPP data in Section 5.1.                  */
/================================================*/

libname rawdata "xxx"; /* Replace xxx with the data source directory. */
proc iml;
start invg(t);
return(1/(1+exp(-t)));
finish;
/* Define the likelihood function assuming perfect tests. */
start lkh(betas) global(groupno, groupy, xs, zs, eps, uneqgs, block5,
block7, block9, subjno, indi);
beta0=betas[1]; beta1=betas[2:3];
py0=1-invg(beta0+xs*beta1);
if indi=1 then logpz0=log(eps+py0);
if indi=0 then do;
if uneqgs=0 then logpz0=groupy*log(eps+py0);
if uneqgs=1 then do;
logpz0=(groupy*log(eps+py0[1:14430]))//
(log(eps+py0[14431:subjno])[+]);
end;
end;
pz1=1-exp(logpz0);
return((zs‘*log(eps+pz1)+(1-zs)‘*logpz0)/groupno);
finish;
/* Define the likelihood function accounting for misclassfication. */
start lkhm(betas) global(groupno, groupy, xs, zs, eta, delta, eps, uneqgs,
block5, block7, block9, subjno, indi);
beta0=betas[1]; beta1=betas[2:3];
py0=1-invg(beta0+xs*beta1);
if indi=1 then pz0=py0;
if indi=0 then do;
if uneqgs=0 then pz0=exp(groupy*log(eps+py0));
if uneqgs=1 then do;
pz0=exp((groupy*log(eps+py0[1:14430]))//
(log(eps+py0[14431:subjno])[+]));
end;
end;
pzm1=eta-(eta+delta-1)#pz0;
return((zs‘*log(pzm1)+(1-zs)‘*log(1-pzm1))/groupno);
finish;
/* Read in the data. */
use rawdata.chlam2008v2; read all into hiv; close rawdata.chlam2008v2;
hiv=hiv[loc(hiv[,1]=1), 2:4];
subjno=nrow(hiv);
print subjno;
xall=hiv[,1:2];
yall=hiv[,3];
/* Provide the sensitivity and specificity. */
eta=0.942; delta=0.976;
/* Form grouped responses: there are 962 pools each of size 15, */
/* and one pool contains 11 subjects. */
groupy15=I((subjno-11)/15)@j(1, 15, 1);
groupno15=963;
/* Form grouped responses for random pooling. */
yrall=yall; xrall=xall;
zr1=1-(groupy15*yall[1:14430]=0);
zr2=1-(yall[14431:subjno][+]=0);
zr15all=zr1//zr2;
/* Form grouped responses for homogeneous pooling, first by age then by race. */
both=yall||xall;
call sort(both, {2 3});
yall=both[,1]; xall=both[,2:3];
zh1=1-(groupy15*yall[1:14430]=0);
zh2=1-(yall[14431:subjno][+]=0);
zh15all=zh1//zh2;
parmno=3; eps=1e-30; optn={1 0}; trueb={-1, 0, 0};

/*~~~~~~~~~~~~~~~~~~~~ Compute naive mles ~~~~~~~~~~~~~~~~~~~~*/
/* individual testing */
zs=yall; xs=xall; groupno=subjno; uneqgs=0; indi=1;
call nlpqn(rc, bni, "lkh", trueb, optn);
call nlpfdd(fv, grd, hes, "lkh", bni);
ainvni=ginv(hes);

/* group testing with random pooling */
zs=zr15all; xs=xrall; groupno=groupno15; groupy=groupy15; uneqgs=1; indi=0;
call nlpqn(rc, bnr15, "lkh", trueb, optn);
call nlpfdd(fv, grd, hes, "lkh", bnr15);
ainvnr15=ginv(hes);
/* group testing with homogeneous pooling */
zs=zh15all; xs=xall; groupno=groupno15; groupy=groupy15; uneqgs=1; indi=0;
call nlpqn(rc, bnh15, "lkh", trueb, optn);
call nlpfdd(fv, grd, hes, "lkh", bnh15);
ainvnh15=ginv(hes);

/*~~~~~~~~~~~~~~~~~~~~ Compute non-naive mles ~~~~~~~~~~~~~~~~~~~~*/
/* individual testing */
zs=yall; xs=xall; groupno=subjno; uneqgs=0; indi=1;
call nlpqn(rc, bti, "lkhm", trueb, optn);
call nlpfdd(fv, grd, hes, "lkhm", bti);
ainvti=ginv(hes);
/* group testing with random pooling */
zs=zr15all; xs=xrall; groupno=groupno15; groupy=groupy15; uneqgs=1; indi=0;
call nlpqn(rc, btr15, "lkhm", trueb, optn);
call nlpfdd(fv, grd, hes, "lkhm", btr15);
ainvtr15=ginv(hes);
/* group testing with homogeneous pooling */
zs=zh15all; xs=xall; groupno=groupno15; groupy=groupy15; uneqgs=1; indi=0;
call nlpqn(rc, bth15, "lkhm", trueb, optn);
call nlpfdd(fv, grd, hes, "lkhm", bth15);
ainvth15=ginv(hes);
/* Compute the gradient matrix for individual testing mles. */
grdni=j(subjno, parmno, 0); grdti=grdni;
groupno=subjno; groupy=1; uneqgs=0; indi=1;
do i=1 to subjno;
zs=yall[i]; xs=xall[i,];
call nlpfdd(fv, grd, useless, "lkh", bni);
grdni[i, ]=grd;
call nlpfdd(fv, grd, useless, "lkhm", bti);
grdti[i, ]=grd;
end;

