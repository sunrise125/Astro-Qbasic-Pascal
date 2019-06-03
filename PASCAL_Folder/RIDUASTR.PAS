Program RIDUASTR;    { RIDUZIONE ASTROMETRICA di un oggetto celeste
                     Tipologia del file dati (nome.DAT) -> Es. A.DAT
GSC2305:786                  1 33  7.420    35 40 55.700   2.895167  4.146875
GSC2305:1126                 1 33  0.780    35 41 26.700   1.894045  4.582708
GSC2305:1074                 1 33  6.060    35 40  9.100   2.666091  3.563272
GSC2305:572                  1 33  0.430    35 37 50.000   1.737662  1.850245
GSC2301:401                  1 32 49.920    35 36 57.500   0.127000  1.247400
GSC2301:693                  1 33 11.120    35 35 57.800   3.319234  0.364718
GSC2301:359                  1 33 22.800    35 37 22.100   5.148568  1.358823
(3908) 1980 PA                                             5.267800  4.028042

--------  lasciare la riga di sopra SEMPRE vuota ---------------------
    Da qui in gi— -> Commenti senza alcuna limitazione
-----------------------------------------------------------------------------
Reference Stars                R.A.            Decl.           Measurement
  Designation                hh mm ss.sss    gg pp ss.sss    x [mm]    y [mm]
-----------------------------------------------------------------------------
              ..... fine file (nome.DAT)

     Note: a) Ampiezza MAX del nome delle stelle e del corpo = 24 caratteri;
              per i dati numerici non Š necessario rispettare la suddetta
              spaziatura, tranne la 44^ col. che d… il segno di ë e risolve
              il caso di ë=-00 gradi
           b) File risultati: (nome.RIS)
}
  Uses Crt;

   Const
         np=20;  rad=pi/180; EPS=1e-10;
   Type
        MatrINTER = Array [1..np] of Integer;
         MatrREAL = Array [1..np] of Real;
        MatrREAL2 = Array [1..np,1..5] of Real;
        MatrREAL3 = Array [1..np,1..np] of Real;
      MatrSTRINGA = Array [1..np] of String[24];
   Var
      fl,ft   : Text;
      w,as    : String;
      obj     : String[24];
      til     : Array [1..np] of String[1];
      x,y,Alfa,Delt,Xg,Yg,ra,de,
      dra,dde,delta             : MatrREAL;
      A                         : MatrREAL2;
      B                         : MatrREAL3;
      nme                       : MatrSTRINGA;
      a0,xp,yp,Xxp,Yyp,
      aw,aa,bb,c,dd,e,f,
      sar,sde,rap,dep           : Real;
      i,nref,band               : Integer;
      s               : Array [1..4] of Real;

      Function ArcoTANG(num,den:Real): Real;
    Var atan: Real;
   Begin
     atan:=arctan(num/den);
      if den<0 then atan:=atan+pi;
      if atan<0 then atan:=atan+2*pi;
     ArcoTANG:=atan;
  End;

      Function ArcoSIN(valore:Real): Real;
   Begin
     ArcoSIN:=arctan(valore/SQRT(1-valore*valore));
  End;

 Function Sw(x:Real): Integer;
  Begin
    if x < 0 then Sw := -1 else Sw := 1;
  End;

 Procedure Stringa;
  Var
      lw: Integer;
      kw,bw,jw,hw            : Real;
      ksw,jsw,hsw               : String;
  Begin
    lw:=Sw(aw);
    aw:=abs(aw);
    kw:=trunc(aw);
    bw:=frac(aw)*60;
    jw:=trunc(bw);
    str(kw:2:0,ksw);
    str(jw:2:0,jsw);
      if band=1 then
         begin
        hw:=trunc(frac(bw)*6000+0.5)/100;
        str(hw:5:2,hsw);
        as:=ksw+' '+jsw+' '+hsw+'';
        end
        else begin
        hw:=trunc(frac(bw)*600+0.5)/10;
        str(hw:4:1,hsw);
        as:=ksw+' '+jsw+' '+hsw+' ';
    if lw<0 then as:='-'+as else as:='+'+as;
      end;

  End;

  Procedure ColoriSchermo;
   Begin
     textcolor(14);textbackground(6);clrscr;
  End;

  Procedure Titolo;
   Begin;
     textcolor(10);
     gotoxy(13,2); writeln('RIDUZIONE ASTROMETRICA - Metodo di Montenbruck-Pfleger');
     gotoxy(13,3); writeln('                                                      ');
     textcolor(14);
  End;

  Procedure SLICE;
   Var
      riga                : String;
      ARs,DEs             : MatrREAL;
      ARh,ARm,DEg,DEp     : MatrINTER;
      segno               : Integer;

    Begin
     gotoxy(24,4); write('Nome del File-Dati (senza ext.)= '); read (w);
     ASSIGN(fl,w+'.dat');
          Reset(fl);
     ASSIGN(ft,w+'.RIS');
          Rewrite(ft);
writeln('     ---------------------------------------------------------------------');
writeln('     Reference Stars            R.A.         Decl.         Measurement   ');
writeln('       Designation           hh mm ss.ss  gg pp ss.s     x [mm]    y [mm]');
writeln('     ---------------------------------------------------------------------');

writeln(ft,'                  A S T R O M E T R I C    R E D U C T I O N              ');
writeln(ft,'     ---------------------------------------------------------------------');
writeln(ft,'     Reference Stars            R.A.         Decl.         Measurement   ');
writeln(ft,'       Designation           hh mm ss.ss  gg pp ss.s     x [mm]    y [mm]');
writeln(ft,'     ---------------------------------------------------------------------');

    nref:=0;
    repeat
       readln(fl,riga);
       nref:=nref+1;        { contatore delle righe UTILI del Dbase }
      til[nref]:= copy(riga,44,1);
      writeln('     '+copy(riga,1,24)+copy(riga,29,11)+' '+copy(riga,44,11)+'   '+copy(riga,60,80));
   writeln(ft,'     '+copy(riga,1,24)+copy(riga,29,11)+' '+copy(riga,44,11)+'   '+copy(riga,60,80));
    until riga='';
       nref:=nref-2;
writeln(ft,'     -------------------------- R I S U L T A T I ------------------------');
          Reset(fl);
    for i:=1 to nref do begin
       readln(fl,nme[i],ARh[i],ARm[i],ARs[i],DEg[i],DEp[i],DEs[i],x[i],y[i]);
       Alfa[i]:=(ARs[i]/3600+ARm[i]/60+ARh[i])*15;
       if  DEg[i]  < 0 then segno:=-1 else segno:=1;
       Delt[i]:=segno*(DEs[i]/3600+DEp[i]/60+ABS(DEg[i]));
             if  (DEg[i]=0) and (til[i]='-') then begin
             Delt[i]:=-Delt[i];
             end;
     end;
      readln(fl,obj,xp,yp);
  Close(fl);

    for i:=1 to nref do begin      {... evita l'Overflow nella Proc. Sistema}
     if x[i]=0 then x[i]:=1e-6; if y[i]=0 then y[i]:=1e-6;
     end;

      readln;  readln;
   End;

   Procedure MATRICE;
   Var i,j: Integer;
   Begin
     for i:=1 to np do begin
      for j:=1 to np do begin
        B[i,j]:=A[i,j];
      end;
     end;
   End;

   Procedure SISTEMA;
   Var i,j,k: Integer;
       p,q,h: Real;
   Begin
     for j:=1 to 3 do begin
      for i:=j+1 to nref do begin
          if B[i,j] <> 0 then begin
               if abs(B[j,j])< EPS*abs(B[i,j]) then begin
                  p:=0; q:=1;
                  B[j,j]:= -B[i,j];
                  B[i,j]:=0;
               end
            else begin
               h:=SQRT(SQR(B[j,j])+SQR(B[i,j]));
               if B[i,j] < 0 then h:=-h;
               p:=B[j,j]/h;
               q:=-B[i,j]/h;
               B[j,j]:=h;
               B[i,j]:=0;
               end;
          end;
         {Calcolo del resto della linea}
          for k:=j+1 to 4 do begin
            h:=p*B[j,k]-q*B[i,k];
            B[i,k]:=q*B[j,k]+p*B[i,k];
            B[j,k]:=h;
          end;
      end;
     end;

         {....sostituzioni}
          for i:=3 downto 1 do begin
           h:=B[i,4];
             for k:=i+1 to 3 do begin
               h:=h+B[i,k]*s[k];
             end;
           s[i]:=-h/B[i,i];
          end;

   End;


   Procedure CALCOLO;
   Var d0,ar0,de0,c,foclen,scale: Real;
   Begin
     Clrscr;
    a0:=0;  d0:=0;
    ar0:=0;  de0:=0;     {AR,DE Centro lastra}
     for i:=1 to nref do begin
      a0:=a0+alfa[i];
      d0:=d0+delt[i];
       alfa[i]:=alfa[i]*rad;
       delt[i]:=delt[i]*rad;
     end;

      a0:=a0/nref; d0:=d0/nref;
      ar0:=a0/15;  de0:=d0;
      a0:=a0*rad; d0:=d0*rad;  { à,ë Centro Lastra in RAD}

   { Calcolo delle Coordinate Standard X[n], Y[n] }

     for i:=1 to nref do begin
      c:=COS(d0)*COS(delt[i])*COS(alfa[i]-a0)+SIN(d0)*SIN(delt[i]);
      Xg[i]:=-(COS(delt[i])*SIN(alfa[i]-a0))/c;
      Yg[i]:=-(SIN(d0)*COS(delt[i])*COS(alfa[i]-a0)-COS(d0)*SIN(delt[i]))/c;
        A[i,1]:=x[i];
        A[i,2]:=y[i];
        A[i,3]:=1;
        A[i,4]:=-Xg[i];
        A[i,5]:=-Yg[i];
     end;

    MATRICE;
    SISTEMA;
    aa:=s[1];
    bb:=s[2];
     c:=s[3];

     for i:=1 to nref do begin
        A[i,4]:= A[i,5];
     end;

    MATRICE;
    SISTEMA;
    dd:=s[1];
     e:=s[2];
     f:=s[3];

 writeln('                               PLATE  Constants     ');
  writeln('          a= ',aa:14:10,'    b= ',bb:14:10,'    c= ',c:14:10);
  writeln('          d= ',dd:14:10,'    e= ',e:14:10,'    f= ',f:14:10);
 writeln(ft,'                               PLATE  Constants     ');
  writeln(ft,'          a= ',aa:14:10,'    b= ',bb:14:10,'    c= ',c:14:10);
  writeln(ft,'          d= ',dd:14:10,'    e= ',e:14:10,'    f= ',f:14:10);

  {Calcolo delle Coordinate equatoriali delle stelle e STIMA degli errori
   ======================================================================}
writeln('  ----------------------------------------------------------------------------');
writeln('    Designation              R.A.          Decl.     -- RESIDUALS in arcsec --');
writeln('                          hh mm ss.ss   gg pp ss.s     R.A.     Decl.    Tot. ');
writeln('  ----------------------------------------------------------------------------');
writeln(ft,'  ----------------------------------------------------------------------------');
writeln(ft,'    Designation              R.A.          Decl.     -- RESIDUALS in arcsec --');
writeln(ft,'                          hh mm ss.ss   gg pp ss.s     R.A.     Decl.    Tot. ');
writeln(ft,'  ----------------------------------------------------------------------------');

band:=1; aw:=ar0; Stringa;  write('  Plate Center            ',as);
                            write(ft,'  Plate Center            ',as);
band:=2; aw:=de0;      Stringa;  writeln('  ',as);
                                 writeln(ft,'  ',as);
writeln;
writeln(ft);

   sar:=0;  sde:=0;

     for i:= 1 to nref do begin
       Xg[i]:=aa*x[i]+bb*y[i]+c;
       Yg[i]:=dd*x[i]+e*y[i]+f;

       ra[i]:=a0+ArcoTANG(-Xg[i],(COS(d0)-Yg[i]*SIN(d0)));
       de[i]:=ArcoSIN((SIN(d0)+Yg[i]*COS(d0))/SQRT(1+SQR(Xg[i])+SQR(Yg[i])));

       ra[i]:=ra[i]/(rad*15);    { in Ore }
            if ra[i]>24 then ra[i]:=ra[i]-24;
       de[i]:=de[i]/rad;         { in Gradi }


       if nref > 3 then begin
            dra[i]:=3600*(15*ra[i]-alfa[i]/rad)*COS(delt[i]);
            dde[i]:=3600*(de[i]-delt[i]/rad);
            sar:=sar+SQR(dra[i]);
            sde:=sde+SQR(dde[i]);
            delta[i]:=SQRT(SQR(dra[i]) + SQR(dde[i]));
            dra[i]:=dra[i];
            end
       else begin dra[i]:=0; dde[i]:=0; delta[i]:=0; end;

band:=1; aw:=ra[i]; Stringa;  write('  '+nme[i]+as);
                           write(ft,'  '+nme[i]+as);

band:=2; aw:=de[i];      Stringa;    write('  ',as);
                                  write(ft,'  ',as);
      writeln('  ',dra[i]:7:3,dde[i]:9:2,delta[i]:9:3);
   writeln(ft,'  ',dra[i]:7:3,dde[i]:9:2,delta[i]:9:3);
     end;

     {posizione dell'oggetto}
     Xxp:=aa*xp+bb*yp+c;
     Yyp:=dd*xp+e*yp+f;
       rap:=a0+ArcoTANG(-Xxp,(COS(d0)-Yyp*SIN(d0)));
       dep:=ArcoSIN((SIN(d0)+Yyp*COS(d0))/SQRT(1+SQR(Xxp)+SQR(Yyp)));

       rap:=rap/(rad*15);    { in Ore }
            if rap>24 then rap:=rap-24;
       dep:=dep/rad;         { in Gradi }
writeln;
writeln(ft);

       if nref > 3 then begin
         sar:=SQRT(sar/(nref-3));
         sde:=SQRT(sde/(nref-3));
            end
       else begin sar:=0; sde:=0; end;

band:=1; aw:=rap; Stringa;  write('  '+obj+as);
                         ;  write(ft,'  '+obj+as);
band:=2; aw:=dep;      Stringa;  write('  ',as);
      writeln('  ',sar:6:2,sde:10:2);
                                 write(ft,'  ',as);
      writeln(ft,'  ',sar:6:2,sde:10:2);
writeln('  ----------------------------------------------------------------------------');
writeln(ft,'  ----------------------------------------------------------------------------');
     foclen:=1/SQRT(abs(aa*e-dd*bb));
     scale:=206264.8/foclen;
  writeln('  Focal Length [mm]:       ',foclen:7:2,'          Image Scale ["/mm]:  ',scale:6:2);
writeln('  ----------------------------------------------------------------------------');
  writeln(ft,'  Focal Length [mm]:       ',foclen:7:2,'          Image Scale ["/mm]:  ',scale:6:2);
writeln(ft,'  ----------------------------------------------------------------------------');

    Close(ft);
   End;

{-------------------------------------------------------------------------------
   INIZIO Programma}
   Begin
       ColoriSchermo;
       Titolo;
       SLICE;
       CALCOLO;
     End.
   {FINE Programma}

