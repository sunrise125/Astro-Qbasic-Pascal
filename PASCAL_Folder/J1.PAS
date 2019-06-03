Program KEYDON;   {TRASFERIMENTO ELLITTICO ELIOCENTRICO INTERPLANETARIO
temp. J.PAS        ‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹

          MEMO: caricare nella Directory di lavoro:
                1) Il file di interfaccia scheda grafica EGAVGA.BGI
                2) I files FONT di caratteri *.CHR
        (altrimenti la parte grafica del programma non va in esecuzione)
}
  Uses  Crt,Graph;
   Const  rad=pi/180; kgauss=0.01720209895; uv0=29.78469169;
          ua=149597870; muS=1;
          np=8; {n.ro Pianeti}

  Type  MatrREAL = Array[1..np] of Real;
        MatrSTRI = Array[1..np] of String[8];

   Var
        x,y,xc,yc, x1,y1,x2,y2, stile0123,spessore01,colore : Word;
        chp: Char;
        Rpia,mu,semmag,massrec: MatrREAL;
                   Nomi,LimInf: MatrSTRI;

        Rpa,Rar,mupa,muar,rinflpa,rinflar,semmagpa,semmagar,
        massrecpa,massrecar,quotapa,
        a,e,p,r1,r2, ka, kSOI, kaINI, cost, t,
        Vp,Va,f1,f2,r1seg,r2seg,v1seg,v2seg,
        deltaV1, ro, deltaV2, alfa,alfaa,
        R0, e1,p1,a1,e2,a2,p2,beta,teta,gamma,gammaa,PerRiv,
        VelioB,VelioBx,VelioBy,Vp1,Vp1x,Vp1y,
        pp1,Vh,Vhx,Vhy,Rparcp1,Vfuga,Vparcp1,
        VelioD,VelioDx,VelioDy,Vp2,Vp2x,Vp2y,
        pp2,Vk,Vkx,Vky,h2,Vparcp2,
        V0, dVhx,dVhy,dVh, xB,yB,xD,yD,xO,yO,fi,fii,psi,psii,f,ff,
        R02,V02    : Real;

{--------------------------- FUNZIONI ----------------------------}
      Function ArcoCOS(x:Real): Real;
   Begin
     ArcoCOS:=pi/2-arctan(x/SQRT(1-x*x));
  End;

      Function ArcoTANG(num,den:Real): Real;
    Var atan: Real;
   Begin
     atan:=arctan(num/den);
      if den<0 then atan:=atan+pi;
      if atan<0 then atan:=atan+2*pi;
     ArcoTANG:=atan;
  End;
{--------------------------- SUB-ROUTINES -------------------------------}
  Procedure ColoriSchermo;
   Begin
     textcolor(14);textbackground(1);clrscr;
  End;

  Procedure Titolo;
   Begin;
     textcolor(11);
     gotoxy(15,2); write('TRASFERIMENTO ELLITTICO ELIOCENTRICO INTERPLANETARIO');
     gotoxy(15,3); write('----------------------------------------------------');
     textcolor(12);
     gotoxy(15,4); write('                               di Veronica Matarazzo');
     textcolor(14);
     gotoxy(15,11); write('------------------------------------------------------------');
     gotoxy(15,20); write('------------------------------------------------------------');
     gotoxy(15,21); write('              Raggi       Ê         a      Masse_Rec.  ÎV1  ');
     textcolor(10);textbackground(0);
     gotoxy(16,23); write('   ATTENZIONE: Venere e Urano hanno rotazione RETROGRADA   ');
     gotoxy(16,24); write('                178¯     98¯  <- Inclin.asse sull''orbita   ');
     textcolor(14);textbackground(1);
  End;

   Procedure LetturaDatiGeoDinamici;
  Var i: Word;
   Begin
     textcolor(14);textbackground(5);
   Nomi[1]:='Mercurio'; Nomi[2]:='Venere  '; Nomi[3]:='Terra   ';
   Nomi[4]:='Marte   '; Nomi[5]:='Giove   '; Nomi[6]:='Saturno ';
   Nomi[7]:='Urano   '; Nomi[8]:='Nettuno ';

   Rpia[1]:=2422; Rpia[2]:=6114; Rpia[3]:=6378; Rpia[4]:=3407;
   Rpia[5]:=71400; Rpia[6]:=57540; Rpia[7]:=24071; Rpia[8]:=24750;

   mu[1]:=21684; mu[2]:=324850; mu[3]:=398600; mu[4]:=42977;
   mu[5]:=126684553; mu[6]:=37911978; mu[7]:=5786122; mu[8]:=6832062;

   semmag[1]:=0.387098; semmag[2]:=0.723332; semmag[3]:=1.000000;
   semmag[4]:=1.523688; semmag[5]:=5.202561; semmag[6]:=9.554747;
   semmag[7]:=19.218140; semmag[8]:=30.109570;

   massrec[1]:=6023600; massrec[2]:=408523.7; massrec[3]:=332946.0;
   massrec[4]:=3098708; massrec[5]:=1047.3; massrec[6]:=3497.9;
   massrec[7]:=22903; massrec[8]:=19412;

   LimInf[1]:=' >1.20  '; LimInf[2]:=' >2.98  '; LimInf[3]:=' >3.25  ';
   LimInf[4]:=' >1.44  '; LimInf[5]:=' >17.4  '; LimInf[6]:=' >10.6  ';
   LimInf[7]:=' >6.40  '; LimInf[8]:='  ----  ';

     for i:=1 to np do begin
    gotoxy(16,11+i);write(i:2,'-',Nomi[i],Rpia[i]:7:0,mu[i]:11:0,semmag[i]:11:6,massrec[i]:11:1,LimInf[i]);
   end;
     textcolor(14);textbackground(1);
  End;

   Procedure InputDati;
 Var nro1,nro2: Integer;
   Begin

repeat
gotoxy(15,6); write('Scelta del pianeta di partenza [1..8] -> ');readln(nro1);
 Rpa:= Rpia[nro1]; mupa:=mu[nro1]; semmagpa:=semmag[nro1];
 massrecpa:=massrec[nro1];
 until nro1 in [1..8];

repeat
gotoxy(15,7); write('Scelta del pianeta di arrivo   [1..8] -> ');readln(nro2);
 Rar:= Rpia[nro2]; muar:=mu[nro2]; semmagar:=semmag[nro2];
 massrecar:=massrec[nro2];
 until (nro2 in [1..8]) and (nro2>nro1);

{ gotoxy(15,8); write('Quota orb. di parcheggio del 1¯ pianeta (in Km)= ');read(h1);}
 gotoxy(15,9);write('Impulso ÎV1 [km/s]= '); readln(DeltaV1);
 gotoxy(45,9);write('Impulso ÎV2 [km/s]= '); readln(DeltaV2);
  clrscr;
  End;

   Procedure SferaInfl;
 Begin
         {Calcolo raggi sfere d'influenza dei pianeti selezionati}
  rinflpa:=ua*semmagpa/exp(0.4*ln(massrecpa));     { in Km }
  rinflar:=ua*semmagar/exp(0.4*ln(massrecar));     { in Km }
 End;

   Procedure ParametriEllisse;   { Sistema eliocentrico }
 Begin
   r1:= semmagpa; r2:= semmagar;           {in UA}
   quotapa:=1/10*rinflpa;
   Rparcp1:=Rpa+quotapa;             {raggio dell'orbita di parcheggio}
   a:=(r1+r2)/2;
   e:=(r2-r1)/(r1+r2);
   p:=a*(1-sqr(e));

   Vp:=sqrt(muS*((2/r1)-(1/a)))*uv0;       {in Km/s}
   Va:=sqrt(muS*((2/r2)-(1/a)))*uv0;       {in Km/s}

   f1:=rinflpa/(ua*r1);     {anomalia del punto B in rad}
   f2:=pi+rinflar/(ua*r2);  {anomalia del punto D in rad}

   r1seg:=p/(1+e*COS(f1));   {rB in UA}
   r2seg:=p/(1+e*COS(f2));   {rD in UA}

   v1seg:=sqrt(muS*((2/r1seg)-(1/a)))*uv0;   {vB in Km/s}
   v2seg:=sqrt(muS*((2/r2seg)-(1/a)))*uv0;   {vD in Km/s}

   PerRiv:=2*pi*(1/kgauss)*sqrt(a*sqr(a)/muS);{Periodo di rivoluzione in giorni}

 End;


   Procedure IperbolePartenza;  { Sistema planetocentrico }
 Var  aux,aux2:Real;
 Begin
{ Componenti della velocitÖ eliocentrica della sonda del punto B }
  VelioBx:=-sqrt(muS/p)*sin(f1)*uv0;
  VelioBy:= sqrt(muS/p)*(e+cos(f1))*uv0;
  VelioB:=sqrt(SQR(VelioBx)+SQR(VelioBy));

  pp1:=r1;              {parametro conica pianeta1= raggiop1}
{ Componenti della velocitÖ eliocentrica del pianeta di partenza }
  Vp1x:=-sqrt(muS/pp1)*sin(f1)*uv0;
  Vp1y:= sqrt(muS/pp1)*cos(f1)*uv0;
  Vp1:=sqrt(SQR(Vp1x)+SQR(Vp1y));

  Vparcp1:=sqrt(mupa/Rparcp1);
  Vfuga:=sqrt(2*mupa/Rparcp1);
  V0:=Vparcp1+DeltaV1;


  Vh:=sqrt(V0*V0-Vfuga*Vfuga);

  e1:=sqr(V0/Vparcp1)-1;    {eccentricitÖ dell'iperbole}
  p1:=rparcp1*(1+e1);            {      p      dell'iperbole}
  a1:=p1/(SQR(e1)-1);            {      a      dell'iperbole}
  aux:=sqrt((a1*a1*(e1*e1-1))/(rinflpa*(2*a1+rinflpa)));
  fi:=pi/2+arcoCOS(aux);
  aux2:=V0/Vparcp1;
  psi:=arctan(aux2*sqrt(sqr(aux2)-2));
  f:=fi-psi;
      {Coordinate di B nel sistema planetocentrico}
  xB:=(r1seg*COS(f1)-r1)*ua;   {in km}
  yB:= r1seg*SIN(f1)*ua;       {in km}
  gamma:=ArcoTANG(yB,xB)-f;

  teta:=pi-psi+gamma;         { angolo di orientamento Vh}
  Vhx:=Vh*COS(teta);
  Vhy:=Vh*SIN(teta);

  {Componenti della Correzione ÎVB}
  dVhx:=VelioBx-Vp1x-Vhx;
  dVhy:=VelioBy-Vp1y-Vhy;
  dVh:=sqrt(SQR(dVhx)+SQR(dVhy));

     xO:=Rparcp1*COS(gamma);   {in km}
     yO:=Rparcp1*SIN(gamma);   {in km}
end;



   Procedure IperboleArrivo;  { Sistema planetocentrico }
 Var  aux,aux2:Real;
  Begin
{ Componenti della velocitÖ eliocentrica della sonda del punto D }
  VelioDx:=-sqrt(muS/p)*sin(f2)*uv0;
  VelioDy:= sqrt(muS/p)*(e+cos(f2))*uv0;
  VelioD:=sqrt(SQR(VelioDx)+SQR(VelioDy));

  pp2:=r2;              {parametro conica pianeta2= raggiop2}
{ Componenti della velocitÖ eliocentrica del pianeta di arrivo }
  Vp2x:=-sqrt(muS/pp2)*sin(f2)*uv0;
  Vp2y:= sqrt(muS/pp2)*cos(f2)*uv0;
  Vp2:=sqrt(SQR(Vp2x)+SQR(Vp2y));

  Vkx:=VelioDx-Vp2x;  {veloc.relativa sonda-pian 2 al lim della sfera infl}
  Vky:=VelioDy-Vp2y;
  Vk:=sqrt(SQR(Vkx)+SQR(Vky));

  if DeltaV2 < Vk/sqrt(2) then begin
     textcolor(142);textbackground(4);
   gotoXY(15,12);write(' E R R O R E: E'' impossibile la cattura al rendez-vous ');
   gotoXY(15,13);write('                   Aumentare la frenata ÎV2= ',DeltaV2:4:2,' km/s ');
     textcolor(14);textbackground(1);
   chp:=readkey;
   halt;
   end;

  Vparcp2:=DeltaV2+sqrt(2*SQR(DeltaV2)-SQR(Vk));
  R02:=muar/SQR(Vparcp2);   {raggio dell'orbita di parcheggio}

  e2:=SQR((Vparcp2+DeltaV2)/Vparcp2)-1;
  a2:=R02/(e2-1);
  p2:=a2*(SQR(e2)-1);
      {Coordinate di D nel sistema planetocentrico}
  xD:=(r2seg*COS(f2)+r2)*ua;   {in km}
  yD:= r2seg*SIN(f2)*ua;       {in km}


  if R02<=Rar then
  Begin
     textcolor(142);textbackground(4);
gotoXY(20,12);write(' E R R O R E: La sonda colpisce il pianeta!      ');
gotoXY(20,13);write('              Diminuire la frenata ÎV2=',DeltaV2:4:2,' km/s ');
     textcolor(14);textbackground(1);
   gotoXY(14,8);write('All''arrivo il raggio dell''orbita di parcheggio ä R=',R02:5:0,' Km');
   gotoXY(14,9);write('che ä minore del raggio del pianeta di arrivo  Rp2=',Rar:5:0,' Km');
   chp:=readkey;
     halt;
  end;

  V02:=Vparcp2+DeltaV2;
  aux:=sqrt((a2*a2*(e2*e2-1))/(rinflar*(2*a2+rinflar)));
  fii:=pi/2+arcoCOS(aux);
  aux2:=V02/Vparcp2;
  psii:=arctan(aux2*sqrt(sqr(aux2)-2));
  ff:=fii-psii;
  gammaa:=ArcoTANG(yD,xD)+ff;

 end;

    Procedure Risultati;
   Begin
     textcolor(14);textbackground(5);
     gotoxy(32,1); writeln('R I S U L T A T I');
  gotoxy(2,3); writeln('Partenza');
     textcolor(14);textbackground(1);
     gotoxy(13,3);write('Angolo Fase ‚= ',gamma/rad:9:4,'¯');
  gotoxy(11,4);
  writeln('  V0=',V0:6:3,'  -> Vparc=',Vparcp1:7:3,' + ÎV1=',DeltaV1:7:3,' km/s');
  gotoxy(11,5);
  writeln('  VB=',VelioB:6:3,'   (  VBx=',VelioBx:7:3,'   VBy=',VelioBy:7:3,' ) km/s');
  gotoxy(11,6);
  writeln(' Vp1=',Vp1:6:3,'   ( Vp1x=',Vp1x:7:3,'  Vp1y=',Vp1y:7:3,' ) km/s');
  gotoxy(11,8);
  writeln('  Vh=',Vh:6:3,'   (  Vhx=',Vhx:7:3,'   Vhy=',Vhy:7:3,' ) km/s');
     textcolor(14);textbackground(4);
  gotoxy(63,9);writeln('Correz. ÎVh');
     textcolor(14);textbackground(1);
  gotoxy(11,9);
  writeln(' ÎVh=',dVh:6:3,'   ( ÎVhx=',dVhx:7:3,'  ÎVhy=',dVhy:7:3,' ) km/s');

  gotoxy(11,11);
  writeln(' Elem_Orbitali Iperb. part.: e=',e1:9:7,'  p=',p1:7:0,' km');
  gotoxy(53,12);
  writeln('a=',a1:7:0,' km');
  gotoxy(11,12);
  writeln(' xB=',xB:7:0,' km     yB=',yB:7:0,' km');
{--------------------------------------------------------}
     textcolor(14);textbackground(5);
  gotoxy(2,14); writeln('Trasfer.');
     textcolor(14);textbackground(1);
  gotoxy(12,14); writeln('Ellisse:  a=',a:8:6,' UA   e=',e:9:7,'   p=',p:8:6,' UA');
  gotoxy(12,15); writeln('        T/2=',PerRiv/2:8:3,' giorni');
{--------------------------------------------------------}
     textcolor(14);textbackground(5);
  gotoxy(2,16); writeln(' Arrivo ');
     textcolor(14);textbackground(1);
  gotoxy(11,16);
  writeln('  VD=',VelioD:6:3,'   (  VDx=',VelioDx:7:3,'   VDy=',VelioDy:7:3,' ) km/s     ÎV2=',DeltaV2:6:3);
  gotoxy(11,17);
  writeln(' Vp2=',Vp2:6:3,'   ( Vp2x=',Vp2x:7:3,'  Vp2y=',Vp2y:7:3,' ) km/s');
  gotoxy(11,18);
  writeln('  Vk=',Vk:6:3,'   (  Vkx=',Vkx:7:3,'   Vky=',Vky:7:3,' ) km/s');

     textcolor(14);textbackground(4);
     gotoxy(63,18);writeln('ÎV2>',Vk/sqrt(2):4:2,' [=Vk/˚2]');
     textcolor(14);textbackground(1);

  gotoxy(11,20);
  writeln('  Vparc=',Vparcp2:7:3,' km/s     R0=',R02:5:0,' km    Quota=',R02-rar:5:0,' km');
  gotoxy(11,22);
  writeln(' Elem_Orbitali Iperb.arrivo: e=',e2:9:7,'  p=',p2:7:0,' km');
  gotoxy(53,23);
  writeln('a=',a2:7:0,' km');
  gotoxy(11,23);
  writeln(' xD=',xD:7:0,' km     yD=',yD:7:0,' km');
  chp:=readkey;  {Blocca lo schermo}
 end;

  Procedure GRAFICOiperbolePartenza;
 Label vai;
   Var
      ch,qq   : String;
      gd,gm,Spost,spost1,spost2,kk,Nparti,FattoreIngr : Integer;
   Begin
    gd:=detect;initgraph(gd,gm,'');
    if graphresult <> grOk then halt(1);
     Spost:=0;   Spost1:=0;   Spost2:=0;
     FattoreIngr:=8;
     kaINI:=60;   {Ingrandimento INIZIALE}
     ka:=kaINI;
     kSOI:=300;Nparti:=1000;
vai:
    cost:=kSOI*ka/ua;
    SetFillStyle(1,1);      { Campitura [0..11] e Colore SFONDO }
    FloodFill(0,0,12);
    SetColor(14);
    SetTextStyle(1,0,1);
    OutTextXY(380,5,'IPERBOLE DI PARTENZA');
    xc:=GetmaxX div 2+Spost; yc:=GetmaxY div 2+Spost2;

     SetFillStyle(1,12);
     SetTextStyle(0,0,1); {Font,Direz.,Size }
    SetLineStyle(0,0,1);
    Rectangle(xc-315,yc+190,xc-75,yc+233);
      SetColor(10);
      OutTextXY(xc-310,yc+195,'+ => Ingrandisce d=> Destra  ');
      OutTextXY(xc-310,yc+208,'- => Rimpicciol. s=> Sinistra');
      OutTextXY(xc-310,yc+220,'USCITA=> Frecce  o=> RESET   ');
      SetColor(14);
    SetLineStyle(3,0,1);
    line(-2*xc-spost1,yc,2*xc+spost1,yc);
  {-----------------------------------------------------------------}
      SetColor(11);
     SetTextStyle(0,0,1); {Font,Direz.,Size}

 FillEllipse(xc,yc,round(Rpa*cost),round(Rpa*cost));{DisegnoPianPart}
 FillEllipse(xc+round(xB*cost),yc-round(yB*cost),3,3);  {Punto B}
    Circle(xc,yc,round(rinflpa*cost));   { Sfera d'influenza}

      SetColor(10);
 Circle(xc,yc,round(Rparcp1*cost));  {Orbita parcheggio}

      SetColor(12);
    SetLineStyle(0,0,1);

 line(xc,yc,xc+round(rparcp1*COS(gamma)*cost),yc-round(rparcp1*SIN(gamma)*cost));

      SetColor(14);
                                 {Iperbole di Partenza}
     alfa:=gamma;
      {1¯ pixel grafico a inizio orbita}
moveto(xc+round(rparcp1*COS(alfa)*cost),yc-round(rparcp1*SIN(alfa)*cost));

     for kk:=0 to (Nparti-1) do begin
        ro:=(p1/(1+e1*COS(alfa-gamma)));
            x:=xc+round(cost*ro*cos(alfa));
            y:=yc-round(cost*ro*sin(alfa));
   Lineto(x,y); {collegam.rettilineo tra le coordinate attuali x,y
                    dello schermo e quelle precedentem.memorizzate}
         Delay(10);
         alfa:=kk/Nparti*f+gamma;
     end;

   Lineto(xc+round(xB*cost),yc-round(yB*cost));

          {----------------------------------}
  repeat
     ch:=readkey;
       if ch='+' then
         begin
         ClearDevice;
          ka:=ka*FattoreIngr;
          goto vai;
          end;

       if ch='-' then
         begin
         ClearDevice;
          ka:=ka/FattoreIngr;
          goto vai;
         end;

       if ch='s' then        {muove figura a sinistra}
         begin
         ClearDevice;
          spost:=spost-100;
          spost1:=spost1+400; {artificio per asse focale}
          goto vai;
          end;

       if ch='d' then        {muove figura a destra}
         begin
         ClearDevice;
          spost:=spost+100;
          goto vai;
          end;

       if ch='a' then       {muove figura in alto}
         begin
         ClearDevice;
          spost2:=spost2-100; FattoreIngr:=4;
          goto vai;
          end;

       if ch='b' then       {muove figura in basso}
         begin
         ClearDevice;
          spost2:=spost2+100; FattoreIngr:=4;
          goto vai;
          end;

       if ch='o' then        {Reset}
         begin
           ka:=kaINI;   {Ingrandimento INIZIALE}
         ClearDevice;
          ka:=kaINI; Spost:=0; Spost1:=0; Spost2:=0;
          goto vai;
          end;
  until keypressed;
    closegraph;
    restoreCrtMode;
  End;


  PROCEDURE PlotOrbita(a,e,teta:real);       {serve per il grafico}
     Begin
           ro:=ka*a*(1-e*e)/(1+e*cos(teta));
            x:=xc+round(ro*cos(teta));
            y:=yc-round(ro*sin(teta));
            Lineto(x,y); {collegam.rettilineo tra le coordinate attuali x,y
                          dello schermo e quelle precedentem.memorizzate}
    end;

   PROCEDURE RaggioVettore(a,e,teta:real;     {serve per il grafico}
         colore,stile0123,spessore13:word);
   Begin
     SetColor(colore);
     SetLineStyle(stile0123,4,spessore13);
     ro:=ka*a*(1-e*e)/(1+e*cos(teta*rad));
     x:=xc+round(ro*cos(teta*rad));
     y:=yc-round(ro*sin(teta*rad));
     Line(xc,yc,x,y);
   end;

  Procedure GRAFICOellisseTrasferimento;
 Label vai;
 Const FattoreIngr=1.4;
   Var
      ch   : String;
      gd,gm,Spost,spost1   : Integer;
      rsole         : Real;
   Begin
    gd:=detect;initgraph(gd,gm,'');
    if graphresult <> grOk then halt(1);
     Spost:=0;
     Spost1:=0;
     ka:=100;   {Ingrandimento INIZIALE}
vai:
    SetFillStyle(1,1);      { Campitura [0..11] e Colore SFONDO }
    FloodFill(0,0,12);
    SetColor(14);
    SetTextStyle(1,0,1);
    OutTextXY(380,5,'ORBITA DI TRASFERIMENTO');
    xc:=GetmaxX div 2+Spost; yc:=GetmaxY div 2;

     SetFillStyle(1,12);
     SetTextStyle(0,0,1); {Font,Direz.,Size }
    SetLineStyle(0,0,1);
    Rectangle(xc-315,yc+190,xc-75,yc+233);
      SetColor(10);
      OutTextXY(xc-310,yc+195,'+ => Ingrandisce d=> Destra  ');
      OutTextXY(xc-310,yc+208,'- => Rimpicciol. s=> Sinistra');
      OutTextXY(xc-310,yc+220,'USCITA=> Frecce  o=> RESET   ');
      SetColor(14);
    rsole:=5;
    SetLineStyle(3,0,1);
    line(-2*xc-spost1,yc,2*xc+spost1,yc);
    FillEllipse(xc,yc,round(rsole),round(rsole));
     SetTextStyle(0,0,1); {Font,Direz.,Size}
     OutTextXY(xc-2,yc+10,'S');
  {-----------------------------------------------------------------}
      SetColor(11);

         Circle(xc,yc,round(r1*ka));    {Disegno Orbita Pianeta di Partenza}
   RaggioVettore (a,e,f1/rad,11,0,1); {a,e,anom.in¯,colore,stile,spessore
                                  stile 0=cont;1=tratt  spess.norm=1,3=triplo}
         Circle(xc+round(r1*ka),yc,round(r1*f1*ka));

      SetColor(12);

         Circle(xc,yc,round(r2*ka));      {Disegno Orbita Pianeta di arrivo}
         RaggioVettore (a,e,f2/rad,12,0,1);
         Circle(xc-round(r2*ka),yc,round(r2*(f2-pi)*ka));

                { Disegno Orbita SONDA per 0 Û teta Û „ step 0.005 rad }
      SetColor(15);
      teta:=0;
      moveto(xc+round(r1*ka),yc); {1¯ pixel grafico a inizio orbita}

     while teta< PI do begin
        PlotOrbita(a,e,teta);
        Delay(10);
      teta:=teta+0.005; {a intervalli di circa 0.3¯}
     end;
          {----------------------------------}
  repeat
     ch:=readkey;
       if ch='+' then
         begin
         ClearDevice;
          ka:=ka*FattoreIngr;
          goto vai;
          end;

       if ch='-' then
         begin
         ClearDevice;
          ka:=ka/FattoreIngr;
          goto vai;
         end;

       if ch='s' then        {muove figura a sinistra}
         begin
         ClearDevice;
          spost:=spost-100;
          spost1:=spost1+400; {artificio per asse focale}
          goto vai;
          end;

       if ch='d' then        {muove figura a destra}
         begin
         ClearDevice;
          spost:=spost+100;
          goto vai;
          end;

       if ch='o' then        {Reset}
         begin
           ka:=100;   {Ingrandimento INIZIALE}
         ClearDevice;
          ka:=150; Spost:=0; Spost1:=0;
          goto vai;
          end;
  until keypressed;
    closegraph;
    restoreCrtMode;
  End;


  Procedure GRAFICOiperboleArrivo;
 Label vai;
   Var
      ch,qq   : String;
      gd,gm,Spost,spost1,spost2,kk,Nparti,FattoreIngr : Integer;
   Begin
    gd:=detect;initgraph(gd,gm,'');
    if graphresult <> grOk then halt(1);
     Spost:=0;   Spost1:=0;   Spost2:=0;
     FattoreIngr:=8;
     kaINI:=60;   {Ingrandimento INIZIALE}
     ka:=kaINI;
     kSOI:=300;Nparti:=1000;
vai:
    cost:=kSOI*ka/ua;
    SetFillStyle(1,1);      { Campitura [0..11] e Colore SFONDO }
    FloodFill(0,0,12);
    SetColor(14);
    SetTextStyle(1,0,1);
    OutTextXY(380,5,'IPERBOLE DI ARRIVO');
    xc:=GetmaxX div 2+Spost; yc:=GetmaxY div 2+Spost2;

     SetFillStyle(1,13);
     SetTextStyle(0,0,1); {Font,Direz.,Size}
    SetLineStyle(0,0,1);

    Rectangle(xc-315,yc+190,xc+40,yc+233);
      SetColor(10);
      OutTextXY(xc-310,yc+195,'+ => Ingrandisce  d=> Destra      a=> Alto ');
      OutTextXY(xc-310,yc+208,'- => Rimpicciol.  s=> Sinistra    b=> Basso ');
      OutTextXY(xc-310,yc+220,'o => RESET      Frecce=> USCITA ');
      SetColor(14);

    SetLineStyle(3,0,1);
    line(-2*xc-spost1,yc,2*xc+spost1,yc);

  {-----------------------------------------------------------------}
      SetColor(11);
     SetTextStyle(0,0,1); {Font,Direz.,Size}

 FillEllipse(xc,yc,round(Rar*cost),round(Rar*cost));{DisegnoPianArrivo}
 FillEllipse(xc+round(xD*cost),yc-round(yD*cost),2,2); {Punto D}
    Circle(xc,yc,round(rinflar*cost));   {Sfera d'influenza}

      SetColor(10);
 Circle(xc,yc,round(R02*cost));  {Orbita parcheggio}

      SetColor(14);
                            {Iperbole di Arrivo}

   SetLineStyle(0,0,1);

     alfaa:=gammaa;

     for kk:=0 downto -Nparti do begin
        ro:=(p2/(1+e2*COS(alfaa-gammaa)));
            x:=xc+round(cost*ro*cos(alfaa));
            y:=yc-round(cost*ro*sin(alfaa));
          Lineto(x,y);
         Delay(10);
         alfaa:=kk/Nparti*ff+gammaa;
     end;

   Lineto(xc+round(xD*cost),yc-round(yD*cost));

         { ---------------------------------- }
  repeat
     ch:=readkey;
       if ch='+' then
         begin
         ClearDevice;
          ka:=ka*FattoreIngr;
          goto vai;
          end;

       if ch='-' then
         begin
         ClearDevice;
          ka:=ka/FattoreIngr;
          goto vai;
         end;

       if ch='s' then        {muove figura a sinistra}
         begin
         ClearDevice;
          spost:=spost-100;
          spost1:=spost1+400; {artificio per asse focale}
          goto vai;
          end;

       if ch='d' then        {muove figura a destra}
         begin
         ClearDevice;
          spost:=spost+100;
          goto vai;
          end;

       if ch='a' then       {muove figura in alto}
         begin
         ClearDevice;
          spost2:=spost2-100; FattoreIngr:=4;
          goto vai;
          end;

       if ch='b' then      {muove figura in basso}
         begin
         ClearDevice;
          spost2:=spost2+100; FattoreIngr:=4;
          goto vai;
          end;

       if ch='o' then        {Reset}
         begin
           ka:=kaINI;  {Ingrandimento INIZIALE}
         ClearDevice;
          ka:=kaINI; Spost:=0; Spost1:=0; Spost2:=0;
          goto vai;
          end;
  until keypressed;
    closegraph;
    restoreCrtMode;
  End;

{---------------------- INIZIO Programma ------------------------------}
   Begin

    ColoriSchermo; Titolo;
    LetturaDatiGeoDinamici;
     InputDati;
     SferaInfl;
      ParametriEllisse;
      IperbolePartenza;
      IperboleArrivo;
     Risultati;
      GRAFICOiperbolePartenza;
      GRAFICOellisseTrasferimento;
      GRAFICOiperboleArrivo;
End.

  {FINE Programma}

