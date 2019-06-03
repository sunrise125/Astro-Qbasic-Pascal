Program EASYHART;  { MOTO PERTURBATO degli N-CORPI risolto l'integratore
                     numerico di Edgar EVERHART (1974) - Universit… di
                     Denver (Colorado)
                                  ------------
    Evoluzione di un corpo celeste a traiettoria ellittica (ASTEROIDE)
    Calcolo automatico dei vettori Pos/Vel del corpo, degli 8 pianeti maggiori
     del sistema solare (Plutone escluso) e della Luna
    Sistema BARICENTRALE - 11 CORPI: (Sole + 8 Pianeti + Luna) + Corpo celeste }

  Uses  Crt;

   Const
       rad=pi/180; kgauss=0.01720209895; ialpha=-1;
       coseps=0.917482137087; sineps=0.397776982902;
       np=11; ZERO=0; HALF=0.5; ONE=1; SR=1.4;

   Label Ripeti;
   Type
        MatrINTER = Array [1..np] of Integer;
         MatrREAL = Array [1..np] of Real;
        MatrREAL2 = Array [1..np*3] of Real;
        MatrREAL3 = Array [1..3] of Real;
        MatrREAL4 = Array [1..21] of Real;
        MatrREAL5 = Array [1..7,1..np*3] of Real;
   Var
      ft,fl                    : Text;
      Sts                      : array[0..13] of Char;
      x0,v0,x,v,f1,fj,y,z      : MatrREAL2;
      am,XS,YS,ZS,VXS,VYS,VZS,AM0       : MatrREAL;
      xx,vv,xt,vt: MatrREAL3;
      c,d,r    : MatrREAL4;
      b,g,e,bd : MatrREAL5;
      w,u      : Array [1..7] of Real;
      nw       : Array [1..8] of Integer;
      h        : Array [1..8] of Real;
      xgeo     : Array [1..3] of Real;
      nomeastero,cata,tm1,tm2    : String;
      il,nv,nb,nclass,LL,ind,
      na,LA,LB,LC,LD,LE,ncount, WTOP, BDTOP, KTOP,
      NI,M,jd,ix,n,k,l,j,i,tt,isom,FLG  : Integer;

     ncl,npq,nes,nsf,nper,npera          : Boolean;
      tf,xl,tm,PW,DIR,tdif,ww,W1,
      TP,SS,T,T2,TVAL,S,Q,A,TEMP,GK,
      HV,rc,rt,dgeolim,
      anom,om,lnodo,incl,eccen,asem, ti       : Real;

  iq,colo,color,xmedq,x1q,y1q,x2q,y2q,uq,xxu,yyu   : Integer;
  ssu,blk,var1,var2 : String;
  scelta  : String[1];
    stri  : String[2];
  col,coli : Array[1..4] of Integer;
  labe : Array[1..4] of String[29];
  ch   : Char;
  lab  : Array[1..4] of Char;
   order : Integer;
 {//////////////////////////////////////////////////////}
       dummy: String;
       gioq,mesq,annq,mesqq,annqq,ap,bp,cp,dp,kk: Real;
       nome: String[20];
       MassaTot,
       Xbar,Ybar,Zbar,VXbar,VYbar,VZbar,
       sx,sy,sz,svx,svy,svz
        : Real;
 {//////////////////////////////////////////////////////}
       giornou,imeseu,iannou,dgeo,rgeo,
          x1,x2,x3,alfa,delt,
          ialfaore,ialfamin,alfasec, idgra,idmin,deltasec,
          alfafraz,deltafraz,Ve,alfat,asemi,sigma,ec,
          whi,whj,whk,wh,wincl,wn,wnodo,ttt,ei,ej,ek,wperi,
          teta,te,Ebig,wanomedia,amotog
                         : Real;
          qqq            : String;
          isegno                : Integer;

      Function ArcoTANG(num,den:Real) : Real;
   Var atan : Real;
   Begin
     atan:=arctan(num/den);
      if den<0 then atan:=atan+pi;
      if atan<0 then atan:=atan+2*pi;
     ArcoTANG:=atan;
  End;

      Function ArcoCOS(valore:Real): Real;
   Begin
     ArcoCOS:=pi/2-arctan(valore/SQRT(ABS(1-valore*valore)));
  End;

      Function ArcoSIN(valore:Real): Real;
   Begin
     ArcoSIN:=arctan(valore/SQRT(ABS(1-valore*valore)));
  End;

   Procedure FORCE(VAR x,F1: MatrREAL2);
   Var
      nm,jm,lm,km              : Integer;
      AAA                      : Real;
      Rh       : Array [1..np,1..np] of Real;
      R1       : MatrREAL;
   Begin
   FOR nm:= 1 TO nb-1 do begin
       jm:= (nm - 1) * 3 + 1;
       na:=nm+1;
     for lm:=na to nb do begin
       km:=(lm-1)*3+1;
  AAA:=SQRT(SQR(x[jm]-x[km]) + SQR(x[jm+1]-x[km+1]) + SQR(x[jm+2]-x[km+2]));
       Rh[nm,lm]:=Exp(3*ln(AAA));
       Rh[nm,lm]:=1/Rh[nm,lm];
       Rh[lm,nm]:=Rh[nm,lm];
      end;
   end;

   for nm:=1 to nb do begin
     jm:=(nm-1)*3+1;
     F1[jm]:=0;  F1[jm+1]:=0;   F1[jm+2]:=0;
    for lm:=1 to nb do begin
      if lm <> nm then begin
       km:=(lm-1)*3+1;
     F1[jm  ]:=F1[jm  ]-AM[lm]*(x[jm  ]-x[km  ])*Rh[nm,lm];
     F1[jm+1]:=F1[jm+1]-AM[lm]*(x[jm+1]-x[km+1])*Rh[nm,lm];
     F1[jm+2]:=F1[jm+2]-AM[lm]*(x[jm+2]-x[km+2])*Rh[nm,lm];
     end;
    end;
   end;
  END; {Fine FORCE}

   Procedure GG2DATA(tm:Real);
   Var
      fu,Zu,alfau,aau,Bu,Cu,Du,eeu       : Real;
   Begin
      Zu:= INT(tm+0.5);
      fu:= tm + 0.5 - Zu;
      alfau:= INT((Zu - 1867216.25) / 36524.25);
      aau:= Zu + 1 + alfau - INT(alfau / 4);
      Bu:= aau + 1524;
      Cu:= INT((Bu - 122.1) / 365.25);
      Du:= INT(365.25 * Cu);
      eeu:= INT((Bu - Du) / 30.6001);
      giornou:= Bu - Du - INT(30.6001 * eeu) + fu;
      IF eeu < 13.5 THEN imeseu:=eeu-1 ELSE imeseu:=eeu-13;
      IF imeseu > 2.5 THEN iannou:=Cu-4716 ELSE iannou:=Cu-4715;
  End;

   Procedure RisultatiFinali(tm,rc,rt:Real;nome,stri:string;
                             NV: Integer;
                             VAR ft:Text; xx,vv,xt,vt: Array of Real);
    Var
          tm1,tm2         :String;

    Begin

     textcolor(12);
         writeln('              RISULTATI FINALI (in Coordinate ECLITTICHE Eliocentriche)');
     textcolor(14);
      writeln(ft,'                                Radau-'+stri);
      writeln(ft,'        RISULTATI FINALI (in Coordinate ECLITTICHE Eliocentriche - J2000.0)');
  write(iannou:6:0,imeseu:3:0,giornou:10:6);
    str(int(tm):7:0,tm1);
    str(frac(tm):7:6,tm2);
    delete(tm2,1,1);
  write('  JD-Finale: ',tm1+tm2);
  textcolor(10); writeln('  Corpo: '+nome); textcolor(14);
  writeln('  Vettore POSIZIONE del Corpo             Vettore VELOCITA` del Corpo');
  writeln(XX[0]:13:9,XX[1]:13:9,XX[2]:13:9,'³',VV[0]:13:9,VV[1]:13:9,VV[2]:13:9);
  writeln('          r= ',rc:13:9);
  write(ft,iannou:6:0,imeseu:3:0,giornou:10:6);
  write(ft,'  JD-Finale: ',tm1+tm2); writeln(ft,'  Corpo: '+nome);
  writeln(ft,'  Vettore POSIZIONE del Corpo             Vettore VELOCITA` del Corpo');
  writeln(ft,XX[0]:13:9,XX[1]:13:9,XX[2]:13:9,'³',VV[0]:13:9,VV[1]:13:9,VV[2]:13:9);
  writeln(ft,'          r= ',rc:13:9);

      Ve:=sqrt((SQR(vv[0])+SQR(vv[1])+SQR(vv[2])))/kgauss;
      alfat:=SQR(Ve)-2/rc;
      asemi:=-1/alfat;
      sigma:=(xx[0]*vv[0]+xx[1]*vv[1]+xx[2]*vv[2])/kgauss;
      ec:=SQRT(1+alfat*(SQR(rc*Ve)-SQR(sigma)));
      whi:=(xx[1]*vv[2]-xx[2]*vv[1])/kgauss;
      whj:=(xx[2]*vv[0]-xx[0]*vv[2])/kgauss;
      whk:=(xx[0]*vv[1]-xx[1]*vv[0])/kgauss;
      wh:=sqrt(SQR(whi)+SQR(whj)+SQR(whk));
      wincl:=ArcoCOS(whk/wh)/rad;
      wn:=sqrt(SQR(whj)+SQR(whi));
      wnodo:=ArcoCOS(-whj/wn)/rad;
      if whi < 0 then wnodo:=360-wnodo;
      ttt:=SQR(Ve)-1/rc;
      ei:=ttt*xx[0]-sigma*vv[0]/kgauss;
      ej:=ttt*xx[1]-sigma*vv[1]/kgauss;
      ek:=ttt*xx[2]-sigma*vv[2]/kgauss;
      wperi:=ArcoCOS((-whj*ei+whi*ej)/(wn*ec))/rad;
      if ek < 0 then wperi:=360-wperi;
      teta:=ArcoCOS((xx[0]*ei+xx[1]*ej+xx[2]*ek)/(rc*ec));
      if sigma < 0 then teta:=2*pi-teta;
      te:=sqrt((1-ec)/(1+ec));
      Ebig:=2*ArcoTANG(te*SIN(teta/2),COS(teta/2));
      wanomedia:=(Ebig-ec*SIN(Ebig))/rad;
      amotog:=180/pi*kgauss/Exp(1.5*ln(asemi));
      writeln ('                           ELEMENTI OSCULATORI');
      writeln ('                    (Eclittica ed Equinozio: J2000.0)');
      writeln('      a = ',asemi:13:8,' UA      e = ',ec:13:8,'    n = ',amotog:12:8,'ø/g');
      writeln('      i = ',wincl:10:5,'ø     ê = ',wnodo:10:5,'ø     w = ',wperi:10:5,'ø');
      writeln('      M = ',wanomedia:10:5,'ø  Epoca= ',iannou:5:0,imeseu:3:0,giornou:9:5);
      writeln;

      writeln (ft,'                           ELEMENTI OSCULATORI');
      writeln (ft,'                    (Eclittica ed Equinozio: J2000.0)');
      writeln(ft,'      a = ',asemi:13:8,' UA      e = ',ec:13:8,'    n = ',amotog:12:8,'ø/g');
      writeln(ft,'      i = ',wincl:10:5,'ø     ê = ',wnodo:10:5,'ø     w = ',wperi:10:5,'ø');
      writeln(ft,'      M = ',wanomedia:10:5,'ø  Epoca= ',iannou:5:0,imeseu:3:0,giornou:9:5);
      writeln (ft);

  writeln('    Vettore POSIZIONE Terra                 Vettore VELOCITA` Terra');
  writeln(XT[0]:13:9,XT[1]:13:9,XT[2]:13:9,'³',VT[0]:13:9,VT[1]:13:9,VT[2]:13:9);
  writeln('          R= ',rt:13:9);
     textcolor(10);
      dgeo:=SQRT(SQR(XX[0]-xt[0])+SQR(XX[1]-xt[1])+SQR(XX[2]-xt[2]));
  writeln('                      Distanza GEOCENTRICA = ',dgeo:10:8);
     textcolor(14);

  writeln(ft,'    Vettore POSIZIONE Terra                 Vettore VELOCITA` Terra');
  writeln(ft,XT[0]:13:9,XT[1]:13:9,XT[2]:13:9,'³',VT[0]:13:9,VT[1]:13:9,VT[2]:13:9);
  writeln(ft,'          R= ',rt:13:9);
  writeln(ft,'                      Distanza GEOCENTRICA = ',dgeo:10:8);
      x1:=xx[0]-xt[0];
      x2:=(xx[1]*coseps-xx[2]*sineps)-(xt[1]*coseps-xt[2]*sineps);
      x3:=(xx[1]*sineps+xx[2]*coseps)-(xt[1]*sineps+xt[2]*coseps);
      rgeo:=sqrt(SQR(x1)+SQR(x2)+SQR(x3));
      alfa:=ArcoTANG(x2,x1);
      if alfa < 0 then alfa:=alfa+2*pi;
      alfa:=alfa*12/pi;
      delt:=ArcoSIN(x3/rgeo)/rad;
         writeln('              à = ',alfa:14:9,' ore;     ë = ',delt:13:8,' gradi');
      writeln(ft,'              à = ',alfa:14:9,' ore;     ë = ',delt:13:8,' gradi');
      ialfaore:= INT(alfa);
      alfafraz:= (alfa - ialfaore) * 60;
      ialfamin:= INT(alfafraz);
      alfasec:= (alfafraz - ialfamin) * 60;
      if delt >= 0 then isegno:=1
      else isegno:=-1;
      idgra:= INT(ABS(delt));
      deltafraz:= (ABS(delt) - idgra) * 60;
      idmin:= INT(deltafraz);
      deltasec:= (deltafraz - idmin) * 60;
      idgra:= idgra * isegno;
      writeln('              à =  ',ialfaore:3:0,ialfamin:3:0,alfasec:7:3,
         '          ë =  ',idgra:3:0,idmin:3:0,deltasec:6:2);
      writeln(ft,'              à =  ',ialfaore:3:0,ialfamin:3:0,alfasec:7:3,
         '          ë =  ',idgra:3:0,idmin:3:0,deltasec:6:2);
   Close(ft);
    end;

   Procedure RA(order: Integer);
   Label 1000,2000,3000,4000,130,150,
         460,490,170,160,180,341,190,360,3,33,333;
   Begin
  Clrscr;
  GG2DATA(tf);
     textcolor(12); gotoxy(2,1);write('RA'+stri); textcolor(10);
 gotoxy(13,1);write(qqq+'          Fine integrazione: ',
                    iannou:5:0,imeseu:3:0,giornou:9:5);
     textcolor(14);
      {Prima chiamata: inizializzazione dell'integratore}
2000:
    nsf:=False;
    nper:=False;
    npera:=False;
    ncl:=nclass=1;
    npq:=nclass<2;
    nes:=LL<0;
       nv:=nb*3;
    if order=11 then begin
  {    Caricamento Vettore NW   }
 nw[1]:= 0; nw[2]:= 0;  nw[3]:= 1;   nw[4]:= 3;
 nw[5]:= 6; nw[6]:= 10;
  {     Caricamento Vettore H (Gauss-Radau Spacings)       }
h[1]:= 0;
h[2]:= 0.09853508579882643;
h[3]:= 0.30453572664636391;
h[4]:= 0.56202518975261386;
h[5]:= 0.80198658212639183;
h[6]:= 0.96019014294853126;
           end
       ELSE   begin
  {    Caricamento Vettore NW   }
 nw[1]:= 0; nw[2]:= 0;  nw[3]:= 1;   nw[4]:= 3;
 nw[5]:= 6; nw[6]:= 10; nw[7]:= 15;  nw[8]:= 21;
  {     Caricamento Vettore H (Gauss-Radau Spacings)       }
h[1]:= 0;
h[2]:= 0.056262560536922146;
h[3]:= 0.180240691736892364;
h[4]:= 0.352624717113169637;
h[5]:= 0.547153626330555383;
h[6]:= 0.734210177215410531;
h[7]:= 0.885320946839095768;
h[8]:= 0.977520613561287501; {' La somma deve essere 3.73333333333333333 (Verif.OK) }
     end;

 tdif:= tf - ti;
  DIR:= tdif / ABS(tdif);
   xl:= DIR * ABS(xl);

{ Calcolo Costanti nei Vettori W(), U() }
  FOR n:= 2 TO WTOP do begin
  ww:= n + n * n;
  IF NCL THEN ww:= n;
  W[n - 1]:= ONE / ww;
  ww:= n;
  u[n - 1]:= ONE / ww;
  End;
{  Azzeramento Vettori BD(), B() }
  FOR k:= 1 TO NV do begin
  IF NCL THEN v[k]:= ZERO;
  FOR l:= 1 TO BDTOP do begin
  bd[l, k]:= ZERO;
   b[l, k]:= ZERO;
  End; End;

{ Calcolo Costanti nei Vettori C(), D(), R() }
  W1:= HALF;
  IF NCL THEN W1:= ONE;
  C[1]:= -h[2]; D[1]:= h[2]; R[1]:= ONE / (h[3] - h[2]);
  LA:= 1; LC:= 1;
   FOR k:= 3 TO KTOP do begin
    LB:= LA; LA:= LC + 1;
    LC:= nw[k + 1];
    C[LA]:= -h[k] * C[LB];
    C[LC]:= C[LA - 1] - h[k];
    D[LA]:= h[2] * D[LB];
    D[LC]:= -C[LC];
    R[LA]:= ONE / (h[k + 1] - h[2]);
    R[LC]:= ONE / (h[k + 1] - h[k]);
        IF k = 3 THEN GOTO 130;
    FOR l:= 4 TO k  do begin
    LD:= LA + l - 3;
    LE:= LB + l - 4;
    C[LD]:= C[LE] - h[k] * C[LE + 1];
    D[LD]:= D[LE] + h[l - 1] * D[LE + 1];
    R[LD]:= ONE / (h[k + 1] - h[l - 1]);
    End;
130: End;
    SS:= Exp(-LL*ln(10));
       TP:= 0.1 * DIR;
       IF xl <> ZERO THEN TP:= xl;
       IF NES THEN TP:= xl;
       IF TP / tdif > HALF THEN TP:= HALF * tdif;
       ncount:= 0;

{*******************************************************************************
}
3000:
      NI:= 6; { Nro delle ITERAZIONI }
  TM:= ti;
          FORCE(x, F1);

4000:
   if nper then npera:=True;

    FOR k:= 1 TO NV  do begin
     if order=11 then begin
G[1, k]:= b[1, k] + D[1] * b[2, k] + D[2] * b[3, k] + D[4] * b[4, k] + D[7] * b[5, k];
G[2, k]:= b[2, k] + D[3] * b[3, k] + D[5] * b[4, k] + D[8] * b[5, k];
G[3, k]:= b[3, k] + D[6] * b[4, k] + D[9] * b[5, k];
G[4, k]:= b[4, k] + D[10]* b[5, k];
G[5, k]:= b[5, k];
           end
        ELSE begin
G[1, k]:= b[1, k] + D[1] * b[2, k] + D[2] * b[3, k] + D[4] * b[4, k] + D[7] * b[5, k] + D[11] * b[6, k] + D[16] * b[7, k];
G[2, k]:= b[2, k] + D[3] * b[3, k] + D[5] * b[4, k] + D[8] * b[5, k] + D[12]* b[6, k] + D[17] * b[7, k];
G[3, k]:= b[3, k] + D[6] * b[4, k] + D[9] * b[5, k] + D[13]* b[6, k] + D[18]* b[7, k];
G[4, k]:= b[4, k] + D[10]* b[5, k] + D[14]* b[6, k] + D[19]* b[7, k];
G[5, k]:= b[5, k] + D[15]* b[6, k] + D[20]* b[7, k];
G[6, k]:= b[6, k] + D[21]* b[7, k];
G[7, k]:= b[7, k];
     end;
    End;

        T:= TP; T2:= T * T;
        IF NCL THEN T2:= T;
        TVAL:= ABS(T);
{                          CICLO PRINCIPALE }

FOR M:= 1 TO NI do begin                 { Loop 490 }
 FOR j:= 2 TO WTOP do begin                 { Loop 470 }
     jd:= j - 1; S:= h[j]; Q:= S;
            IF NCL THEN Q:= ONE;
       FOR k:= 1 TO NV do begin          { Loop 150 }
     if order=11 then begin
    A:=W[3]*b[3,k]+S*(W[4]*b[4,k]+S*(W[5]*b[5,k]));
           end
        ELSE begin
    A:=W[3]*b[3,k]+S*(W[4]*b[4,k]+S*(W[5]*b[5,k]+S*(W[6]*b[6,k]+S*W[7]*b[7,k])));
        end;
 Y[k]:=x[k]+Q*(T*v[k]+T2*S*(F1[k]*W1+S*(W[1]*b[1,k]+S*(W[2]*b[2,k]+S*A))));
       IF NPQ THEN goto 150;
     if order=11 then begin
    A:=u[3]*b[3,k]+S*(u[4]*b[4,k]+S*(u[5]*b[5,k]));
           end
        ELSE begin
    A:=u[3]*b[3,k]+S*(u[4]*b[4,k]+S*(u[5]*b[5,k]+S*(u[6]*b[6,k]+S*u[7]*b[7,k])));
        end;
   Z[k]:=v[k]+S*T*(F1[k]+S*(u[1]*b[1,k]+S*(u[2]*b[2,k]+S*A)));
150:
  End;

      FORCE(Y, FJ);

       FOR k:= 1 TO NV do begin          { Loop 460 }
         TEMP:= G[jd, k];
         GK:= (FJ[k] - F1[k])/S;

   if (j=1) or (j=2) then begin
      G[1, k]:=GK;
          goto 360;
          End;
   if j=3 then begin
      G[2, k]:=(GK-G[1,k])*R[1];
          goto 360;
          End;
   if j=4 then begin
     G[3, k]:=((GK-G[1,k])*R[2]-G[2,k])*R[3];
          goto 360;
          End;
   if j=5 then begin
     G[4, k]:=(((GK-G[1,k])*R[4]-G[2,k])*R[5]-G[3,k])*R[6];
          goto 360;
          End;
   if j=6 then begin
     G[5, k]:=((((GK-G[1,k])*R[7]-G[2,k])*R[8]-G[3,k])*R[9]-G[4,k])*R[10];
          goto 360;
          End;

     if order=15 then begin
   if j=7 then begin
     G[6, k]:=(((((GK-G[1,k])*R[11]-G[2,k])*R[12]-G[3,k])*R[13]-G[4,k])*R[14]-G[5,k])*R[15];
          goto 360;
          End;
   if j=8 then begin
     G[7, k]:=((((((GK-G[1,k])*R[16]-G[2,k])*R[17]-G[3,k])*R[18]-G[4,k])*R[19]-G[5,k])*R[20]-G[6,k])*R[21];
          End;
     end;
360:   TEMP:= G[jd, k] - TEMP;
        b[jd, k]:= b[jd, k] + TEMP;

   if (j=1) or (j=2) then goto 460;
   if j=3 then begin
     b[1, k]:= b[1, k] + C[1] * TEMP;
          goto 460;
          End;
   if j=4 then begin
     b[1, k]:= b[1, k] + C[2] * TEMP;
     b[2, k]:= b[2, k] + C[3] * TEMP;
          goto 460;
          End;
   if j=5 then begin
     b[1, k]:= b[1, k] + C[4] * TEMP;
     b[2, k]:= b[2, k] + C[5] * TEMP;
     b[3, k]:= b[3, k] + C[6] * TEMP;
          goto 460;
          End;
   if j=6 then begin
     b[1, k]:= b[1, k] + C[7] * TEMP;
     b[2, k]:= b[2, k] + C[8] * TEMP;
     b[3, k]:= b[3, k] + C[9] * TEMP;
     b[4, k]:= b[4, k] + C[10] * TEMP;
          goto 460;
          End;
   if j=7 then begin
     b[1, k]:= b[1, k] + C[11] * TEMP;
     b[2, k]:= b[2, k] + C[12] * TEMP;
     b[3, k]:= b[3, k] + C[13] * TEMP;
     b[4, k]:= b[4, k] + C[14] * TEMP;
     b[5, k]:= b[5, k] + C[15] * TEMP;
          goto 460;
          End;

     if order=15 then begin
   if j=8 then begin
     b[1, k]:= b[1, k] + C[16] * TEMP;
     b[2, k]:= b[2, k] + C[17] * TEMP;
     b[3, k]:= b[3, k] + C[18] * TEMP;
     b[4, k]:= b[4, k] + C[19] * TEMP;
     b[5, k]:= b[5, k] + C[20] * TEMP;
     b[6, k]:= b[6, k] + C[21] * TEMP;
          End;
     end;
460: End;
     End;
  IF NES OR (M<NI) THEN goto 490;
  HV:= ZERO;
   FOR k:= 1 TO NV do begin
     if order=11 then begin
   IF ABS(b[5, k]) > HV THEN HV:= ABS(b[5, k]); { Max tra HV e abs(B(5,k)) }
           end
        ELSE begin
   IF ABS(b[7, k]) > HV THEN HV:= ABS(b[7, k]); { Max tra HV e abs(B(7,k)) }
     end;
   End;
     if order=11 then begin
   HV:= HV * W[5] / Exp(5*ln(TVAL));
           end
        ELSE begin
   HV:= HV * W[7] / Exp(7*ln(TVAL));
     end;
 490: End;
 IF nsf THEN goto 170;
 IF NOT NES THEN TP:= Exp(PW*ln(SS/HV)) * DIR;
 IF NES THEN TP:= xl;
 IF NES THEN goto 160;
 IF TP / T > ONE THEN goto 160;
     TP:= 0.8 * TP;
     ncount:= ncount + 1;
 { NCOUNT Š una sicurezza sul numero dei cicli: NON DEVE superare 10,
   altrimenti il Programma SI FERMA }
IF ncount >= 10 THEN begin
writeln('     Attenzione !! ... NCOUNT Š uguale a ', ncount:3);
goto 1000;
end;
160: nsf:= True;
170: FOR k:= 1 TO NV do begin
     if order=11 then begin
     x[k]:=x[k]+v[k]*T+T2*(F1[k]*W1+b[1,k]*W[1]+b[2,k]*W[2]+b[3,k]*W[3]+b[4,k]*W[4]+b[5,k]*W[5]);
           end
        ELSE begin
     x[k]:=x[k]+v[k]*T+T2*(F1[k]*W1+b[1,k]*W[1]+b[2,k]*W[2]+b[3,k]*W[3]+b[4,k]*W[4]+b[5,k]*W[5]+b[6,k]*W[6]+b[7,k]*W[7]);
     end;
     IF NCL THEN goto 180;
     if order=11 then begin
     v[k]:=v[k]+T*(F1[k]+b[1,k]*u[1]+b[2,k]*u[2]+b[3,k]*u[3]+b[4,k]*u[4]+b[5,k]*u[5]);
           end
        ELSE begin
     v[k]:=v[k]+T*(F1[k]+b[1,k]*u[1]+b[2,k]*u[2]+b[3,k]*u[3]+b[4,k]*u[4]+b[5,k]*u[5]+b[6,k]*u[6]+b[7,k]*u[7]);
     end;
180: End;
{ Mio inserimento: visualizza ad ogni ciclo x(),v(),TM=tempo progressivo }
  gotoxy(8,2);
  writeln('Vettore POSIZIONE                      Vettore VELOCITA`');
  writeln(X[NV-2]-x[1]:13:9,X[NV-1]-x[2]:13:9,X[NV]-x[3]:13:9,'³',
          V[NV-2]-v[1]:13:9,V[NV-1]-v[2]:13:9,V[NV]-v[3]:13:9);

    isom:=8;
   FOR ix:= 1 TO 3 do begin
   xx[ix]:= x[NV - 3 + ix] - x[ix];       { Ritorno al Sistema Eliocentrico }
   vv[ix]:= v[NV - 3 + ix] - v[ix];
   xt[ix]:= x[NV-3*isom + ix] - x[ix];
   vt[ix]:= v[NV-3*isom + ix] - v[ix];
   end;
  TM:= TM + T;
    dgeo:=SQRT(SQR(XX[1]-xt[1])+SQR(XX[2]-xt[2])+SQR(XX[3]-xt[3]));
    str(int(tm):7:0,tm1);
    str(frac(tm):7:6,tm2);
    delete(tm2,1,1);

  writeln('         JD=',tm1+tm2,'  Interv.: ',T:12:7,'    DIST_geo: ',dgeo:13:8);
  GG2DATA(tm);
  writeln(iannou:13:0,imeseu:3:0,giornou:10:6);
      if dgeo < dgeolim then begin
      writeln(ft,'         ',tm1+tm2,'   ',iannou:6:0,imeseu:3:0,giornou:10:6,
                 '    DIST_geo: ',dgeo:13:8);
      end;

      {  Torna se ha finito  }

 IF nper then GOTO 1000;

    FORCE(x, F1);

  IF NES THEN GOTO 341;
  TP:= DIR * Exp(PW*ln(SS/HV));
  IF TP / T > SR THEN TP:= T * SR;

341:  IF NES THEN TP:= xl;
   IF DIR * (TM + TP) < DIR * tf - 0.00000001 THEN goto 190;
   TP:= tf - TM;
    nper:= True;

190: Q:= TP / T;
 FOR k:= 1 TO NV do begin  {Loop 220}
      FOR j:= 1 TO BDTOP do begin
      bd[j, k]:= b[j, k] - E[j, k];
      End;

     if order=11 then begin
   E[1,k]:= Q *     (b[1, k]+2*b[2, k]+3 *b[3, k]+ 4*b[4, k]+ 5*b[5, k]);
   E[2,k]:= Q * Q * (b[2, k]+3*b[3, k]+6 *b[4, k]+10*b[5, k]);
   E[3,k]:= Exp(3*ln(Q))*(b[3, k]+4*b[4, k]+10*b[5, k]);
   E[4,k]:= Exp(4*ln(Q))*(b[4, k]+5*b[5, k]);
   E[5,k]:= Exp(5*ln(Q))*b[5, k];
           end
        ELSE begin
   E[1,k]:= Q *     (b[1, k]+2*b[2, k]+3 *b[3, k]+ 4*b[4, k]+ 5*b[5, k]+ 6*b[6, k]+7*b[7, k]);
   E[2,k]:= Q * Q * (b[2, k]+3*b[3, k]+6 *b[4, k]+10*b[5, k]+15*b[6, k]+21*b[7, k]);
   E[3,k]:= Exp(3*ln(Q))*(b[3, k]+4*b[4, k]+10*b[5, k]+20*b[6, k]+35*b[7, k]);
   E[4,k]:= Exp(4*ln(Q))*(b[4, k]+5*b[5, k]+15*b[6, k]+35*b[7, k]);
   E[5,k]:= Exp(5*ln(Q))*(b[5, k]+6*b[6, k]+21*b[7, k]);
   E[6,k]:= Exp(6*ln(Q))*(b[6, k]+7*b[7, k]);
   E[7,k]:= Exp(7*ln(Q))* b[7, k];
     end;

      FOR l:= 1 TO BDTOP do begin b[l, k]:= E[l, k] + bd[l, k]; End;
 End;
   NI:= 2; { Solo 2 iterazioni per tutte le sequenze dopo la prima }

      goto 4000;
1000:
{          ULTIMA USCITA
           -------------            }

  GG2DATA(tm);

    isom:=8;
   FOR ix:= 1 TO 3 do begin
   xx[ix]:= x[NV - 3 + ix] - x[ix];       { Ritorno al Sistema Eliocentrico }
   vv[ix]:= v[NV - 3 + ix] - v[ix];
   xt[ix]:= x[NV-3*isom + ix] - x[ix];
   vt[ix]:= v[NV-3*isom + ix] - v[ix];
   end;
      rc:=SQRT(SQR(xx[1])+SQR(xx[2])+SQR(xx[3]));
      rt:=SQRT(SQR(xt[1])+SQR(xt[2])+SQR(xt[3]));
   RisultatiFinali(tm,rc,rt,nome,stri,NV,ft,xx,vv,xt,vt);

  End;
   {FINE Procedura RA(order) }


 Procedure MessLine;
  Begin
   gotoxy(1,24);
   TextBackground(7);
   TextColor(4);
   writeln(blk);
  End;

 Procedure Retina;
  Begin
    TextBackGround(6);
    ClrScr;
    ssu:='';
    for uq:=1 to 79 do
     begin
      ssu:=ssu+chr(176);
     end;
    for uq:=1 to 25 do
     begin
      writeln(ssu);
     end;
  End;

 Procedure BoxB;
  Begin
   ssu:='';
   xmedq:=(x1q+x2q) div 2;
   xxu:=x2q-x1q;
   yyu:=y2q-y1q;

   for uq:=1 to xxu-2 do
    begin
     ssu:=ssu+chr(196);
    end;
    TextColor(4);
    TextBackGround(7);

    gotoxy(x1q,y1q);
    writeln(chr(218),ssu,chr(191));
    for uq:=1 to yyu-1 do
     begin
      gotoxy(x1q,wherey);
      writeln(chr(179),chr(179):(xxu-1));
     end;
     gotoxy(x1q,wherey);
     writeln(chr(192),ssu,chr(217));
 End;

 Procedure BordEx;
  begin
   ssu:='';
   xxu:=x2q-x1q;
   yyu:=y2q-y1q;

   for uq:=1 to xxu do
    begin
    ssu:=ssu+chr(176);
    end;
    gotoxy(x1q+1,y2q+1);
    TextBackGround(16);
    Textcolor(8);
    writeln(ssu);
    gotoxy(x2q,y1q+1);
    for uq:=1 to yyu do
     begin
      writeln(chr(176));
      gotoxy(x2q,wherey);
     end;
  End;

 Procedure Piccola;
  begin
   TextBackground(colo);
   TextColor(color);
   gotoxy(27,8+2*(iq-1));
   writeln(labe[iq]);
  end;

 Procedure MenuIn;
  Begin
   lab[1]:='1';
   lab[2]:='2';
   lab[3]:='3';
   labe[1]:='Integ. di EVERHART 11ø ordine';
   labe[2]:='Integ. di EVERHART 15ø ordine';
   labe[3]:='    >>>> U S C I T A  <<<<   ';
   for iq:=1 to 3 do
    begin
     col[iq]:=7;
     coli[iq]:=4;
    end;

   col[1]:=9;
   coli[1]:=14;

   TextColor(9);
   Retina;
   MessLine;
   gotoxy(20,24);
   TextColor(1);
   writeln('Muovere le frecce o selezionare i numeri');
   ssu:='';
   x1q:=23;
   y1q:=6;
   x2q:=57;
   y2q:=14;
   BoxB;
   BordEx;
   textcolor(4);
   textbackground(7);
   gotoxy(xmedq-13,y1q);
   writeln('  Moto PERTURBATO N-Corpi  ');
   gotoxy(xmedq-12,y2q+3);
   writeln(' Ing. Giuseppe MATARAZZO ');
   gotoxy(xmedq-16,y2q-2);
 { textcolor(1);
   writeln('                                ');
   gotoxy(xmedq-16,y2q-1);
   writeln('                                ');
   textcolor(4);                                   }

    for iq:=1 to 3 do
     begin
      colo:=col[iq];
      color:=coli[iq];
      textbackground(7);
      textcolor(0);
      gotoxy(25,8+2*(iq-1));
      writeln(lab[iq]+'.');
      Piccola;
     end;

     gotoxy(75,23);
     iq:=1;

     repeat
     ch:=readkey;
     if ch=#72 then
      begin
       colo:=7;
       color:=4;
       Piccola;
       iq:=iq-1;
       if iq<1 then iq:=3;
       colo:=9;
       color:=14;
       Piccola;
      end;
     if ch=#80 then
      begin
       colo:=7;
       color:=4;
       Piccola;
       iq:=iq+1;
       if iq>3 then iq:=1;
       colo:=9;
       color:=14;
       Piccola;
      end;

      gotoxy(75,23);

     until (ch=#13) or (ch=lab[1]) or (ch=lab[2])
                    or (ch=lab[3]);

     if ch=#13 then str(iq:1,scelta)
     else
     begin
      for iq:=1 to 3 do
       begin
        if ch=lab[iq] then str(iq:1,scelta);
       end;
     end;
  End;

  Procedure ColoriSchermo;
   Begin
     textcolor(14);textbackground(6);clrscr;
  End;

  Procedure Icona;
   Begin                                                                                  ;
writeln;                                                                                  ;
writeln('  ²²²²²²      ²²     ²²²²²²   ²²  ²²   ²²    ²²     ²²     ²²²²²²²²   ²²²²²²²² ');
writeln('  ²²        ²²  ²²   ²²       ²²  ²²   ²²    ²²   ²²  ²²   ²²    ²²      ²²    ');
writeln('  ²²       ²²    ²²  ²²         ²²     ²²    ²²  ²²    ²²  ²²    ²²      ²²    ');
writeln('  ²²²²²    ²²²²²²²²    ²²       ²²     ²²²²²²²²  ²²²²²²²²  ²²²²²²²²      ²²    ');
writeln('  ²²       ²²    ²²      ²²     ²²     ²²    ²²  ²²    ²²  ²²  ²²        ²²    ');
writeln('  ²²       ²²    ²²      ²²     ²²     ²²    ²²  ²²    ²²  ²²   ²²       ²²    ');
writeln('  ²²       ²²    ²²      ²²     ²²     ²²    ²²  ²²    ²²  ²²    ²²      ²²    ');
writeln('  ²²²²²²   ²²    ²²  ²²²²²²     ²²     ²²    ²²  ²²    ²²  ²²     ²²     ²²    ');
writeln;
writeln(' ÚÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄ¿');
writeln(' ³ MOTO PERTURBATO degli N-CORPI risolto con l''integratore di Edgar EVERHART  ³');
writeln(' ³ - Universit… di Denver (Colorado-USA) - Evoluzione di un corpo celeste a   ³');
writeln(' ³ traiettoria ellittica (ASTEROIDE). Vettori Posizione e Velocit… del Corpo, ³');
writeln(' ³ degli 8 pianeti (no Plut.) e della Luna all''inizio integr. ricavati con il ³');
writeln(' ³ pacchetto applicativo CERES dell''ITA di S.Pietroburgo.                     ³');
writeln(' ³                                Set.99                                    ³');
writeln(' ³                                                                            ³');
writeln(' ³      Input Dati:  FILE ESTERNO: EASYHART.DAT (Vettori Pos/Vel dei corpi)  ³');
writeln(' ³                                                                            ³');
writeln(' ³Output Risultati:  SU VIDEO: Pos/Vel e dist.geocentr. corpo ad ogni step di³');
writeln(' ³                    integraz. e finali. Elementi Orbitali di fine integraz. ³');
writeln(' ³                   FILE ESTERNO:(Nome).RIS; Dist.geoc.parziali + Elem.Orb. ³');
writeln(' ÀÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÙ');
readln;
  End;

   Procedure DataFineIntegrazione(var jdtop:Real);
 Var m: WORD;
     y: INTEGER;
  g,mm,yy,a,b,kk,c,d: REAL;
  Begin
gotoxy(17,4); write('ÚÄÄÄÄÄÄÄÄÄ Data di FINE Integrazione ÄÄÄÄÄÄ¿');
gotoxy(17,5); write('³                                          ³');
gotoxy(17,6); write('³                                          ³');
gotoxy(17,7); write('³                                          ³');
gotoxy(17,8); write('³                                          ³');
gotoxy(17,9); write('³                                          ³');
gotoxy(17,10);write('ÀÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÙ');
repeat gotoxy(25,5);write('  Giorno.dec = ');readln(g);
      until (g<32) and (g>=0);
repeat gotoxy(25,6);write('        Mese = ');readln(m);
      until m in [1..12];
       gotoxy(25,7);write('        Anno = ');readln(y);

       yy:=y; mm:=m;
          IF m <= 2 then
           begin
           y:=y-1; m:=m+12;
           end;
              b:=0; kk:=yy+mm/100+g/10000;
          IF kk > 1582.101499 then    {Data >= 15.10.1582--->Rif.Gregoriana}
           begin
          a:=int(y/100); b:=2 - a + int(a/4);
           end;
              c:=int(365.25*y);
              d:=int(30.6001*(m+1));
      jdtop:=b+c+d+g+1720994.5;
 End;

   Procedure PassoDistGeoLim(var xxl,ddgeolim:Real);
   Begin
gotoxy(17,4); write('ÚÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄ¿');
gotoxy(17,5); write('³                                          ³');
gotoxy(17,6); write('³                                          ³');
gotoxy(17,7); write('³                                          ³');
gotoxy(17,8); write('³                                          ³');
gotoxy(17,9); write('³                                          ³');
gotoxy(17,10);write('ÀÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÙ');
   gotoxy(22,6);write('       (0 = Int. AUTOMATICO)     ');
   gotoxy(22,5);write('Interv. Integrazione [Giorni] = ');readln(xxl);
   gotoxy(22,8);write(' Limite Dist.Geocentrica [UA] = ');readln(ddgeolim);
   End;

   Procedure FileEsterno;
 Label Skip1,Skip2;
   Begin

   FLG:=FLG+1;
            if FLG>1 then goto Skip1;

   DataFineIntegrazione(tf);
   PassoDistGeoLim(xl,dgeolim);
      clrscr;
      writeln('     Giorno Giuliano di Inizio Integrazione  : ',ti:15:5);
      writeln('     Giorno Giuliano di   Fine Integrazione  : ',tf:15:5);
      writeln('     Intervallo di integraz. (0=automatico)  : ',xl:12:2);
      writeln('     LIMITE della distanza geocentrica       : ',dgeolim:12:2);
Skip1:

        if FLG>1 then goto Skip2;
      writeln;
      writeln('          Vettori POSIZIONE |-- Coord.BARICENTR.--| Vettori VELOCITA`');
     for i:=1 to np do
      Begin
      writeln(XS[i]:13:9,YS[i]:13:9,ZS[i]:13:9,'³',VXS[i]:13:9,VYS[i]:13:9,VZS[i]:13:9);
     End;
      readln;
Skip2:

  {.......... trasporto dati nei Vettori (x0) e (v0) ............}
     for i:=1 to np do
      Begin
      j:=(i-1)*3+3;
      x0[j-2]:=XS[i]; v0[j-2]:=VXS[i];
      x0[j-1]:=YS[i]; v0[j-1]:=VYS[i];
      x0[j  ]:=ZS[i]; v0[j  ]:=VZS[i];
     End;

  { Selezione dei CORPI Attivi }
     il:=0;

  FOR j:=1 to np do begin
      il:= il + 1;
       FOR k:= 1 TO 3 do
        Begin
           tt:=(il - 1) * 3 + k;
        x[(il - 1) * 3 + k]:= x0[(j - 1) * 3 + k];
        v[(il - 1) * 3 + k]:= v0[(j - 1) * 3 + k];
       End;
       am[il]:= am0[j];
     END;

 nb:=il; {Numero corpi attivi}

  { SETTAGGIO dei parametri d'integrazione }
      nv:=nb*3;
      nclass:=-2;
      if ABS(xl)> 1e-10 then LL:=-LL;

     ASSIGN(ft,copy(nome,13,8)+'.R'+stri);
          Rewrite(ft);
     writeln(ft,'     FILTRO su DIST_Geoc minore di: '
                ,dgeolim:5:2,' UA   Interv.[0=Automatico] = ',xl:6:2,' g');

   End; {Fine Proc. FileEsterno}

{-------------------------------------------------------------------------------
   INIZIO Programma}
   Begin
   ColoriSchermo; Icona;



{--------- Passaggio dal Sistema Eclittico Eliocentrico ----------------------
           al Sistema Eclittico Baricentrale ---------------------------------}
  sx:=0; sy:=0;  sz:=0;  svx:=0; svy:=0;  svz:=0;
  MassaTot:=0;

  {Parametri Gravitazionali -> æ_sole= kgaussý= 0.01720209895ý UA^3/gý }
am0[1]:= 2.95912208285591e-4;  { Sole }
        {masse inverse dei pianeti (fonte: DE403-JPL) al denominatore}
am0[2]:= am0[1]/6023600;           am0[3]:= am0[1]/408523.71; { æ: Mer,Ven }
am0[4]:= am0[1]/332946.0486301812; am0[5]:= am0[1]/3098708;   { æ: Ter,Mar }
am0[6]:= am0[1]/1047.3486;         am0[7]:= am0[1]/3497.898;  { æ: Gio,Sat }
am0[8]:= am0[1]/22902.98;          am0[9]:= am0[1]/19412.24;  { æ: Ura,Net }
am0[10]:= am0[1]/27068708.52707218; am0[11]:= am0[1]/1e20;    { æ: Luna, Corpo di massa piccolissima}

   for i:=1 to np do begin
     MassaTot:= MassaTot+am0[i];
   end;

    ASSIGN (fl,'EASYHART.DAT');
    Reset(fl);

   for i:=1 to np-1 do begin
        {Lettura pos/vel. elioc. Pian.perturb.+Corpo}
   readln (fl,nome);
   readln (fl,dummy);
   readln (fl,dummy);
   readln (fl,dummy);
   readln (fl,dummy);
   readln (fl,annq);
   readln (fl,dummy);
   readln (fl,mesq,gioq,XS[i],YS[i],ZS[i]);
   readln (fl,VXS[i],VYS[i],VZS[i]);
   end;
   Close(fl);

    {JD della data di inizio integrazione}
       annqq:=annq; mesqq:=mesq;
          IF mesq <= 2 then begin
           annq:=annq-1; mesq:=mesq+12;
           end;
              bp:=0; kk:=annqq+mesqq/100+gioq/10000;
          IF kk > 1582.101499 then    {Data >= 15.10.1582--->Rif.Gregoriana}
           begin
          ap:=int(annq/100); bp:=2 - ap + int(ap/4);
           end;
              cp:=int(365.25*annq);
              dp:=int(30.6001*(mesq+1));
      ti:=bp+cp+dp+gioq+1720994.5;

{-------------------}

   for i:=np downto 2 do begin
   XS[i]:= XS[i-1];  YS[i]:= YS[i-1];  ZS[i]:= ZS[i-1];
   VXS[i]:= VXS[i-1];  VYS[i]:= VYS[i-1];  VZS[i]:= VZS[i-1];
   end;

            {Pos/Vel Sole, elioc.}
   XS[1]:=0; YS[1]:=0; ZS[1]:=0;   VXS[1]:=0; VYS[1]:=0;  VZS[1]:=0;

          {Calcolo Baricentro Pos/Vel Sistema Solare}
   for i:=1 to np do begin
    sx:=sx+am0[i]*XS[i];  sy:=sy+am0[i]*YS[i];  sz:=sz+am0[i]*ZS[i];
    svx:=svx+am0[i]*VXS[i];  svy:=svy+am0[i]*VYS[i];  svz:=svz+am0[i]*VZS[i];
   end;

       Xbar:=sx/MassaTot;  Ybar:=sy/MassaTot;   Zbar:=sz/MassaTot;
       VXbar:=svx/MassaTot;  VYbar:=svy/MassaTot;   VZbar:=svz/MassaTot;

 {Pos/Vel Sole+Pianeti in coord.eclittiche baricentrali}
  For i:= 1 to np do begin
   XS[i]:=XS[i]-Xbar;  YS[i]:=YS[i]-Ybar;   ZS[i]:=ZS[i]-Zbar;
   VXS[i]:=VXS[i]-VXbar;  VYS[i]:=VYS[i]-VYbar;   VZS[i]:=VZS[i]-VZbar;
{  writeln(i:2,XS[i]:12:8,YS[i]:12:8,ZS[i]:12:8,'³',VXS[i]:12:9,VYS[i]:13:9,VZS[i]:13:9);}
  end;
 {........................ Come UNIQUE99 .............}
   blk:='';
   for uq:=1 to 79 do
    begin
     blk:=blk+chr(32);
    end;

     FLG:=0;

 Ripeti:
   MenuIn;
   ch:=' ';

   if scelta='1' then
    Begin
     textcolor(14);textbackground(6);clrscr;
    order:=11; str(order:2,stri);  LL:=8;   PW:= 1 / 7;
                                 WTOP:=6; BDTOP:=5; KTOP:=5;
    FileEsterno;
    RA(order);
    readln;
    end;

   if scelta='2' then
    Begin
     textcolor(14);textbackground(6);clrscr;
    order:=15; str(order:2,stri);  LL:=10;  PW:= 1 / 9;
                                 WTOP:=8; BDTOP:=7; KTOP:=7;
    FileEsterno;
    RA(order);
    readln;
    end;

   if scelta='3' then
    Begin
    exit;
    end;
 Goto Ripeti;

     End.
   {FINE Programma}
