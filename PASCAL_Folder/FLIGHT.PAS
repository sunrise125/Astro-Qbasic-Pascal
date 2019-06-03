Program Flight;   {TRASFERIMENTO ELLITTICO ELIOCENTRICO INTERPLANETARIO}

  Uses  Crt,Graph;
   Const  rad=pi/180; kgauss=0.01720209895; uv0=29.78469169;
          ua=149597870; muS=1;
          np=8; {n.ro Pianeti}

  Type  MatrREAL = Array[1..np] of Real;
        MatrSTRI = Array[1..np] of String[8];

   Var  x,y,xc,yc, x1,y1,x2,y2, stile0123,spessore01,colore : Word;
        scelta,chp:    Char;
        Rpia,mu,semmag,massrec:    MatrREAL;
        Nomi,LimInf:    MatrSTRI;
        nro1,nro2:   Integer;
        Rpa,Rar,mupa,muar,rinflpa,rinflar,semmagpa,semmagar,
        massrecpa,massrecar,
        a,e,p,r1,r2, ka, kSOI, kaINI, cost, t,
        Vp,Va,f1,f2,r1seg,r2seg,v1seg,v2seg, r0fb,
        deltaV1,deltaV1min, ro, deltaV2, alfa,alfaa,
        R0, e1,p1,a1,e2,a2,p2,beta,teta,gamma,gammaa,PerRiv,
        VelioB,VelioBx,VelioBy,Vp1,Vp1x,Vp1y,
        pp1,Vh,Vhx,Vhy,h1,Rparcp1,Rparcp2,Vfuga,Vparcp1,
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
 PROCEDURE ColoriSchermo;
   Begin
      textcolor(14);textbackground(1);clrscr;
   End;

 PROCEDURE Copertina;
 var gd,gm:integer;
  Begin
   gd:=detect;initgraph(gd,gm,'');
   if graphresult <> grOk then halt(1);
   SetBkColor(1);SetColor(11);
   SetTextStyle(1,0,3);OutTextXY(10,50,'Trasferimento Ellittico Eliocentrico Interplanetario ');
   SetFillStyle(1,11);Bar(10,90,595,90);SetColor(13);
   SetTextSTYLE(1,0,2);OutTextXY(370,130,'di Veronica Matarazzo');SetColor(13);
   SetTextSTYLE(1,0,2);SetColor(15);
   OutTextXY(110,240,' Simulazione del moto di una sonda ');
   OutTextXY(110,270,' da un pianeta interno ad uno pi— ');
   OutTextXY(110,300,'    esterno del Sistema Solare ');
   SetTextStyle(1,0,1);SetColor(14);OutTextXY(330,420,' Premere INVIO per continuare ');readln;
   CloseGraph;
   RestoreCrtMode;
  End;

 PROCEDURE Titolo;
   Begin;
     textcolor(11);
     gotoxy(15,2); write('TRASFERIMENTO ELLITTICO ELIOCENTRICO INTERPLANETARIO');
     gotoxy(15,3); write('ÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄ');
     textcolor(12);
     gotoxy(15,4); write('                               di Veronica Matarazzo');
     textcolor(14);
     gotoxy(15,11); write('ÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄ-');
     gotoxy(15,20); write('ÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄ-');
     gotoxy(15,21); write('              Raggi       æ         a      Masse_Rec. ');
     textcolor(10);textbackground(0);
     gotoxy(16,23); write('   ATTENZIONE: Venere e Urano hanno rotazione RETROGRADA   ');
     gotoxy(16,24); write('                178ø     98ø  <- Inclin.asse sull''orbita   ');
     textcolor(14);textbackground(1);
   End;

 PROCEDURE LetturaDatiGeoDinamici;
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

    for i:=1 to np do begin
     gotoxy(16,11+i);write(i:2,'-',Nomi[i],Rpia[i]:7:0,mu[i]:11:0,semmag[i]:11:6,massrec[i]:11:1,LimInf[i]);
    end;
     textcolor(14);textbackground(1);
  End;

   PROCEDURE InputDati;
   Begin

repeat
gotoxy(15,6); write('Scelta del pianeta di partenza [1..8] -> ');readln(nro1);
 Rpa:= Rpia[nro1]; mupa:=mu[nro1]; semmagpa:=semmag[nro1]; massrecpa:=massrec[nro1];
 until nro1 in [1..8];

repeat
gotoxy(15,7); write('Scelta del pianeta di arrivo   [1..8] -> ');readln(nro2);
 Rar:= Rpia[nro2]; muar:=mu[nro2]; semmagar:=semmag[nro2]; massrecar:=massrec[nro2];
 until (nro2 in [1..8]) and (nro2>nro1);

 gotoxy(15,8); write('Quota orb. di parcheggio del 1ø pianeta (in Km)= ');read(h1);
 deltaV1min:=(sqrt(2)-1)*sqrt(mupa/(Rpa+h1));
 repeat
 gotoxy(15,9);write('Impulso ëV1 [>',deltaV1min:4:3,']= '); readln(deltaV1);
 until deltaV1>deltaV1min;
  clrscr;
  End;

   PROCEDURE SferaInfl;
 Begin
         {Calcolo raggi sfere d'influenza dei pianeti selezionati}
  rinflpa:=ua*semmagpa/exp(0.4*ln(massrecpa));     { in Km }
  rinflar:=ua*semmagar/exp(0.4*ln(massrecar));     { in Km }
 End;

   PROCEDURE ParametriEllisse;   { Sistema eliocentrico }
 Begin
   r1:= semmagpa; r2:= semmagar;           {in UA}
   Rparcp1:=Rpa+h1;             {raggio dell'orbita di parcheggio}
   a:=(r1+r2)/2;
   e:=(r2-r1)/(r1+r2);
   p:=a*(1-sqr(e));

   Vp:=sqrt(muS*((2/r1)-(1/a)))*uv0;       {in Km/s}
   Va:=sqrt(muS*((2/r2)-(1/a)))*uv0;       {in Km/s}

   f1:=rinflpa/(ua*r1);     {anomalia del punto B in rad}
   f2:=pi-rinflar/(ua*r2);  {anomalia del punto D in rad}

   r1seg:=p/(1+e*COS(f1));   {rB in UA}
   r2seg:=p/(1+e*COS(f2));   {rD in UA}

   v1seg:=sqrt(muS*((2/r1seg)-(1/a)))*uv0;   {vB in Km/s}
   v2seg:=sqrt(muS*((2/r2seg)-(1/a)))*uv0; Vky:=-Vky; Vkx:=Vkx+Vky/10;

   PerRiv:=2*pi*(1/kgauss)*sqrt(a*sqr(a)/muS);{Periodo di rivoluzione in giorni}

 End;


   PROCEDURE IperbolePartenza;  { Sistema planetocentrico }
 Var  aux,aux2:Real;
 Begin
{ Componenti della velocit… eliocentrica della sonda del punto B }
  VelioBx:=-sqrt(muS/p)*sin(f1)*uv0;
  VelioBy:= sqrt(muS/p)*(e+cos(f1))*uv0;
  VelioB:=sqrt(SQR(VelioBx)+SQR(VelioBy));

  pp1:=r1;              {parametro conica pianeta1= raggiop1}
{ Componenti della velocit… eliocentrica del pianeta di partenza }
  Vp1x:=-sqrt(muS/pp1)*sin(f1)*uv0;
  Vp1y:= sqrt(muS/pp1)*cos(f1)*uv0;
  Vp1:=sqrt(SQR(Vp1x)+SQR(Vp1y));

  Vparcp1:=sqrt(mupa/Rparcp1);
  Vfuga:=sqrt(2*mupa/Rparcp1);
  V0:=Vparcp1+DeltaV1;


  Vh:=sqrt(V0*V0-Vfuga*Vfuga);

  e1:=sqr(V0/Vparcp1)-1;    {eccentricit… dell'iperbole}
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

  {Componenti della Correzione ëVB}
  dVhx:=VelioBx-Vp1x-Vhx;
  dVhy:=VelioBy-Vp1y-Vhy;
  dVh:=sqrt(SQR(dVhx)+SQR(dVhy));

end;


   PROCEDURE IperboleArrivo;  { Sistema planetocentrico }
 Var  aux,aux2,efb,dfb,Vkxa,Vkya,Vka,prodsca :Real;
  Begin
{ Componenti della velocit… eliocentrica della sonda del punto D }
  VelioDx:=-sqrt(muS/p)*sin(f2)*uv0;
  VelioDy:= sqrt(muS/p)*(e+cos(f2))*uv0;
  VelioD:=sqrt(SQR(VelioDx)+SQR(VelioDy));

  pp2:=r2;              {parametro conica pianeta2= raggiop2}
{ Componenti della velocit… eliocentrica del pianeta di arrivo }
  Vp2x:=-sqrt(muS/pp2)*sin(f2)*uv0;
  Vp2y:= sqrt(muS/pp2)*cos(f2)*uv0;
  Vp2:=sqrt(SQR(Vp2x)+SQR(Vp2y));

  Vkx:=VelioDx-Vp2x;  {veloc.relativa sonda-pian 2 al lim della sfera infl}
  Vky:=VelioDy-Vp2y; ParametriEllisse;
  Vk:=sqrt(SQR(Vkx)+SQR(Vky));

  writeln('  Vk=',Vk:6:3,'  ð (  Vkx=',Vkx:7:3,'   Vky=',Vky:7:3,' ) km/s');


      {Coordinate di D nel sistema planetocentrico}
  xD:=(r2seg*COS(f2)+r2)*ua;   {in km}
  yD:= r2seg*SIN(f2)*ua;       {in km}

  writeln(' xD=',xD:7:0,' km     yD=',yD:7:0,' km');

       {Calcolo del FlyBy}
  Vkxa:=-Vkx;  Vkya:=Vky;  Vka:=Vk;
  prodsca:= (Vkx*Vkxa+Vky*Vkya)/(Vk*Vka);
  dfb:= ArcoCOS(prodsca);
  efb:= 1/SIN(dfb/2);
  r0fb:=muar*(efb-1)/SQR(Vk);
  writeln(' r0fb=',r0fb:7:0,' km    efb=',efb:10:6);

       {Parametri iperbole di arrivo: e,a,p e Impulso di Frenata}
R02:=r0fb/10;
  V02:=sqrt(Vk*Vk+2*muar/R02);
  Vparcp2:=sqrt(muar/R02);
  aux2:=V02/Vparcp2;
  e2:=SQR(aux2)-1;
  a2:=R02/(e2-1);
  p2:=a2*(SQR(e2)-1);

  DeltaV2:=V02-Vparcp2;      {impulso di velocit… all'arrivo: frenata}

  writeln;
  writeln(' e2=',e2:10:6,'   a2=',a2:7:0,'  km    p2=',p2:7:0,' km');
  Rparcp2:=R02; h2:=R02-Rar;
  writeln(' R02=',R02:7:0,'   h2=',h2:7:0,'  km');
  writeln(' DeltaV2=',DeltaV2:7:3,' km/s');
  clrscr;

  aux:=sqrt((a2*a2*(e2*e2-1))/(rinflar*(2*a2+rinflar)));
  fii:=pi/2+arcoCOS(aux);
  aux2:=V02/Vparcp2;
  psii:=arctan(aux2*sqrt(sqr(aux2)-2));
  ff:=fii-psii;
  gammaa:=ArcoTANG(yD,xD)+ff;

 end;

    PROCEDURE Risultati;
   Begin
     textcolor(14);textbackground(5);
     gotoxy(32,1); writeln('R I S U L T A T I');
  gotoxy(2,3); writeln('Partenza');
     textcolor(14);textbackground(1);
     gotoxy(13,3);write('Angolo Fase â= ',gamma/rad:9:4,'ø');
  gotoxy(11,4);
  writeln('  V0=',V0:6:3,'  -> Vparc=',Vparcp1:7:3,' + ëV1=',DeltaV1:7:3,' km/s');
  gotoxy(11,5);
  writeln('  VB=',VelioB:6:3,'  ð (  VBx=',VelioBx:7:3,'   VBy=',VelioBy:7:3,' ) km/s');
  gotoxy(11,6);
  writeln(' Vp1=',Vp1:6:3,'  ð ( Vp1x=',Vp1x:7:3,'  Vp1y=',Vp1y:7:3,' ) km/s');
  gotoxy(11,8);
  writeln('  Vh=',Vh:6:3,'  ð (  Vhx=',Vhx:7:3,'   Vhy=',Vhy:7:3,' ) km/s');
     textcolor(14);textbackground(4);
  gotoxy(63,9);writeln('Correz. ëVh');
     textcolor(14);textbackground(1);
  gotoxy(11,9);
  writeln(' ëVh=',dVh:6:3,'  ð ( ëVhx=',dVhx:7:3,'  ëVhy=',dVhy:7:3,' ) km/s');

  gotoxy(11,11);
  writeln(' Elem_Orbitali Iperb. part.: e=',e1:9:7,'  p=',p1:7:0,' km');
  gotoxy(53,12);
  writeln('a=',a1:7:0,'Km');
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
  gotoxy(2,17); writeln(' Arrivo ');
     textcolor(14);textbackground(1);
  gotoxy(11,17);
  writeln('  VD=',VelioD:6:3,'  ð (  VDx=',VelioDx:7:3,'   VDy=',VelioDy:7:3,' ) km/s');
  gotoxy(11,18);
  writeln('  Vk=',Vk:6:3,'  ð (  Vkx=',Vkx:7:3,'   Vky=',Vky:7:3,' ) km/s  ëV2=',DeltaV2:6:3,' km/s');
  gotoxy(11,20);
  writeln('  Vparc=',Vparcp2:7:3,' km/s   Rparc=',Rparcp2:5:0,' km   Quota=',h2:5:0,' km');
  gotoxy(11,22);
  writeln(' Elem_Orbitali Iperb.arrivo: e=',e2:9:7,'  p=',p2:7:0,' km');
  gotoxy(11,23); writeln(' xD=',xD:7:0,' km     yD=',yD:7:0,' km');
  gotoxy(53,23); writeln('a=',a2:7:0,' km');
  gotoxy(11,24); writeln(' gamma=',gammaa/rad:8:3,'ø');

  chp:=readkey;  {Blocca lo schermo}
 end;

  PROCEDURE GRAFICOiperbolePartenza;
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
    Rectangle(5,430,365,478);
      SetColor(10);
      OutTextXY(15,435,'+ => Ingrandisce  d=> Destra      a=> Alto ');
      OutTextXY(15,450,'- => Rimpicciol.  s=> Sinistra    b=> Basso ');
      OutTextXY(15,465,'o => RESET        Frecce=> CONTINUA ');
      SetColor(14);
    SetLineStyle(3,0,1);
    line(-2*xc-spost1,yc,2*xc+spost1,yc);
  {-----------------------------------------------------------------} 
      SetColor(11);
     SetTextStyle(0,0,1); {Font,Direz.,Size}

 FillEllipse(xc,yc,round(Rpa*cost),round(Rpa*cost));{DisegnoPianPart}
 FillEllipse(xc+round(xB*cost),yc-round(yB*cost),3,3);  {Punto B}
    Circle(xc,yc,round(rinflpa*cost));   {Sfera d'influenza}

      SetColor(10);
 Circle(xc,yc,round(Rparcp1*cost));  {Orbita parcheggio}

      SetColor(5);
    SetLineStyle(0,0,1);

 line(xc,yc,xc+round(rparcp1*COS(gamma)*cost),yc-round(rparcp1*SIN(gamma)*cost));

      SetColor(14);
                                 {Iperbole di Partenza}
     alfa:=gamma;
      {1ø pixel grafico a inizio orbita}
moveto(xc+round(rparcp1*COS(alfa)*cost),yc-round(rparcp1*SIN(alfa)*cost));

     for kk:=0 to (Nparti-1) do begin
        ro:=(p1/(1+e1*COS(alfa-gamma)));
            x:=xc+round(cost*ro*cos(alfa));
            y:=yc-round(cost*ro*sin(alfa));
   Lineto(x,y); {collegam.rettilineo tra le coordinate attuali x,y
                    dello schermo e quelle precedentem.memorizzate}
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

  PROCEDURE GRAFICOellisseTrasferimento;
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

    Rectangle(5,430,300,478);
      SetColor(10);
      OutTextXY(15,435,'+ => Ingrandisce  d=> Destra ');
      OutTextXY(15,450,'- => Rimpicciol.  s=> Sinistra ');
      OutTextXY(15,465,'o => RESET        Frecce=> CONTINUA ');
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
  RaggioVettore (a,e,f1/rad,11,0,1); {a,e,anom.inø,colore,stile,spessore
                                   stile 0=cont;1=tratt  spess.norm=1,3=triplo}
      Circle(xc+round(r1*ka),yc,round(r1*f1*ka));
      SetColor(12);
      Circle(xc,yc,round(r2*ka));      {Disegno Orbita Pianeta di arrivo}
      RaggioVettore (a,e,f2/rad,12,0,1);
      Circle(xc-round(r2*ka),yc,round(r2*(pi-f2)*ka));

                { Disegno Orbita SONDA per 0 ó teta ó ã step 0.005 rad }
      SetColor(15);
      teta:=0;
      moveto(xc+round(r1*ka),yc); {1ø pixel grafico a inizio orbita}

     while teta< PI do begin
        PlotOrbita(a,e,teta);
      teta:=teta+0.005; {a intervalli di circa 0.3ø}
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


  PROCEDURE GRAFICOiperboleArrivo;
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

    Rectangle(5,430,365,478);
      SetColor(10);
      OutTextXY(15,435,'+ => Ingrandisce  d=> Destra      a=> Alto ');
      OutTextXY(15,450,'- => Rimpicciol.  s=> Sinistra    b=> Basso ');
      OutTextXY(15,465,'o => RESET        Frecce=> MENU'' ');
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

      SetColor(14);       {Iperbole di Arrivo}
   SetLineStyle(0,0,1);
     alfaa:=gammaa;

     for kk:=0 downto -Nparti do begin
        ro:=(p2/(1+e2*COS(alfaa-gammaa)));
            x:=xc+round(cost*ro*cos(alfaa));
            y:=yc-round(cost*ro*sin(alfaa));
   Lineto(x,y);
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



PROCEDURE Introduzione;
 var gd,gm:integer;
     ch:char;
  Begin
    gd:=detect;initgraph(gd,gm,'');
    if graphresult <> grOk then halt(1);
    SetGraphMode(3);ClearViewPort;SetViewPort(0,0,GetMaxX,GetMaxY,ClipOn);
    SetColor(13);SetTextStyle(1,0,4);OutTextXY(100,0,' Introduzione al problema ');
    SetBkColor(1);SetFillStyle(1,13);
    Bar(110,40,520,43);
    SetTextStyle(2,0,6);SetColor(11);
    OutTextXY(30,50,'Nel trasferimento di sonde artificiali da un pianeta ad un');
    OutTextXY(30,70,'altro del Sistema Solare Š importante tener presente, tra ');
    OutTextXY(30,90,'i tanti fattori in gioco, quello del consumo energetico.  ');
   OutTextXY(30,110,'                                                          ');
   OutTextXY(30,130,'La traiettoria di minima energia totale si realizza quando');
   OutTextXY(30,150,'il pianeta di partenza si trova al perielio e quello di   ');
   OutTextXY(30,170,'arrivo all''afelio, quando cioŠ l''angolo di trasferimento Š');
   OutTextXY(30,190,'di 180ø. Una situazione ideale come questa Š difficile che');
   OutTextXY(30,210,'si verifichi in realt…, ma alcune deduzioni teoriche che  ');
   OutTextXY(30,230,'essa implica sono di grande ausilio per la comprensione   ');
   OutTextXY(30,250,'dei fenomeni gravitazionali presenti.                     ');
   OutTextXY(30,270,'                                                          ');
   OutTextXY(30,290,'La sonda deve prima sottrarsi all''azione newtoniana della ');
   OutTextXY(30,310,'massa del pianeta di partenza, con un dispendio ëV1 di    ');
   OutTextXY(30,330,'energia prodotta dal combustibile, quindi spostarsi, sotto');
   OutTextXY(30,350,'l''influenza gravitazionale del Sole, verso il pianeta di  ');
   OutTextXY(30,370,'arrivo ed infine penetrare dentro la sua sfera di         ');
   OutTextXY(30,390,'influenza e posizionarsi, consumando in frenata un altro  ');
   OutTextXY(30,410,'ëV2 di carburante, in un''orbita temporanea di parcheggio. ');
   SetTextStyle(1,0,1); SetColor(15);
   OutTextXY(155,450,'Premere INVIO per continuare');
   readln;
   ClearViewPort; SetTextStyle(2,0,6); SetColor(11);
   OutTextXY(30,130,'Scopo di questo lavoro Š quello di simulare una missione  ');
   OutTextXY(30,150,'spaziale come quella descritta considerando, per          ');
   OutTextXY(30,170,'semplificare i calcoli analitici, che le orbite dei       ');
   OutTextXY(30,190,'pianeti siano circolari e complanari e che inoltre il     ');
   OutTextXY(30,210,'trasferimento avvenga da un pianeta interno verso uno     ');
   OutTextXY(30,230,'esterno del Sistema Solare.                               ');
   OutTextXY(30,250,'                                                          ');
   SetTextStyle(1,0,1); SetColor(15);
   OutTextXY(100,450,'Premere una freccia per ritornare al menu''  ');
   ch:=readkey;
    CloseGraph;
    RestoreCrtMode;
 End;

PROCEDURE Automatico;
 var gd,gm:integer;
     ch:char;
  Begin
    gd:=detect;initgraph(gd,gm,'');
    if graphresult <> grOk then halt(1);
    SetGraphMode(3);ClearViewPort;SetViewPort(0,0,GetMaxX,GetMaxY,ClipOn);
    SetColor(15);SetTextStyle(3,0,4);OutTextXY(0,0,' Simulazione mediante dati prestabiliti ');
    SetBkColor(1); Setfillstyle(1,15); Bar(20,50,600,52);

    SetColor(10);SetTextStyle(1,0,2);
    OutTextXY(15,120,' Pianeta di partenza = Terra ');
    OutTextXY(15,170,' Pianeta di arrivo   = Marte ');
    OutTextXY(15,220,' Quota orbita di parcheggio del primo pianeta   = 200 Km ');
    OutTextXY(15,270,' Impulso iniziale DeltaV1 = 3.6 Km/sec ');
    SetColor(15);
    OutTextXY(80,420,' Premere invio per osservare il trasferimento ');
    ch:=readkey;
    CloseGraph;
    RestoreCrtMode;
 {-------------------------------}
    ColoriSchermo;
     LetturaDatiGeoDinamici; clrscr;
    nro1:=3;     nro2:=4;
    h1:=200;     h2:=400;
    DeltaV1:=3.6;
    Rpa:= Rpia[nro1]; mupa:=mu[nro1];
          semmagpa:=semmag[nro1]; massrecpa:=massrec[nro1];
    Rar:= Rpia[nro2]; muar:=mu[nro2];
          semmagar:=semmag[nro2]; massrecar:=massrec[nro2];
    SferaInfl;
     ParametriEllisse;
     IperbolePartenza;
     IperboleArrivo;
     Risultati;
          GRAFICOIperbolePartenza;
          GRAFICOellisseTrasferimento;
          GRAFICOiperboleArrivo;
  End;


  PROCEDURE Avvio;
  Begin
     Titolo;
     LetturaDatiGeoDinamici;
     InputDati;
     SferaInfl;
     ParametriEllisse;
     IperbolePartenza;
     IperboleArrivo;
     Risultati;
         GRAFICOIperbolePartenza;
         GRAFICOellisseTrasferimento;
         GRAFICOiperboleArrivo;
  End;

  PROCEDURE Avvio2;
  Begin
     SferaInfl;
     ParametriEllisse;
     IperbolePartenza;
     IperboleArrivo;
     Risultati;
         GRAFICOIperbolePartenza;
         GRAFICOellisseTrasferimento;
         GRAFICOiperboleArrivo;
  End;

PROCEDURE Generico;
 var gd,gm:integer;
     ch:char;
  Begin
    textcolor(14);textbackground(1);
  gotoxy(15,12); write('-----------------------------------------------------');
  gotoxy(22,13); write('Questi sono i dati reali dei pianeti:');
  gotoxy(15,14); write('-----------------------------------------------------');
  gotoxy(15,15); write(' 1-Mercurio   2422      21684   0.387098  6023600.0  ');
  gotoxy(15,16); write(' 2-Venere     6114     324850   0.723332   408523.7  ');
  gotoxy(15,17); write(' 3-Terra      6378     398600   1.000000   332946.0  ');
  gotoxy(15,18); write(' 4-Marte      3407      42977   1.523688  3098708.0  ');
  gotoxy(15,19); write(' 5-Giove     71400  126684553   5.202561     1047.3  ');
  gotoxy(15,20); write(' 6-Saturno   57540   37911978   9.554747     3497.9  ');
  gotoxy(15,21); write(' 7-Urano     24071    5786122  19.218140    22903.0  ');
  gotoxy(15,22); write(' 8-Nettuno   24750    6832062  30.109570    19412.0  ');
  gotoxy(15,23); write('-----------------------------------------------------');
  gotoxy(15,24); write('              Raggi      æ          a      Masse_Rec.');
textcolor(13);
gotoxy(32,2); write('DATI GENERICI');
textcolor(11);
gotoxy(10,4); write('PIANETA DI PARTENZA                   PIANETA DI ARRIVO');
textcolor(15);
gotoxy(12,5); write('Raggio_1 [Km]= ');readln(Rpa);
gotoxy(50,5); write('Raggio_2 [Km]= ');readln(Rar);
gotoxy(12,6); write('æ1 [Km^3/sý]= ');readln(mupa);
gotoxy(50,6); write('æ2 [Km^3/sý]= ');readln(muar);
gotoxy(7,7); write('Semiasse orbitale_1 [UA]= ');readln(semmagpa);
repeat
gotoxy(45,7); write('Semiasse orbitale_2 [UA]= ');readln(semmagar);
until semmagar>semmagpa;
gotoxy(2,8); write('Massa rec_1(=Msole/Mp1)[adim.]=');readln(massrecpa);
gotoxy(40,9); write('Massa rec_2(=Msole/Mp2)[adim.]=');readln(massrecar);
gotoxy(2,10); write('Quota orbita di parcheggio_1 [Km]= ');readln(h1);
deltaV1min:=(sqrt(2)-1)*sqrt(mupa/(Rpa+h1));
 repeat
gotoxy(7,11);write('Impulso ëV1 [>',deltaV1min:4:3,']= '); readln(deltaV1);
 until deltaV1>deltaV1min;
Avvio2;
  End;


PROCEDURE Simulazione;
 var gd,gm:integer;
     scelta,ch:char;
Begin
  gd:=detect;initgraph(gd,gm,'');
  if graphresult <> grOk then halt(1);
  SetGraphMode(3);ClearViewPort;SetViewPort(0,0,GetMaxX,GetMaxY,ClipOn);
  SetFillStyle(7,7);
  Bar(120,5,500,40); SetColor(15); Rectangle(120,5,500,40); Rectangle(115,0,505,45);
  SetTextStyle(1,0,4); SetColor(14); OutTextXY(125,0,' Simulazione del moto ');
  SetBkColor(9);SetColor(15);
  Setfillstyle(7,12);
  Bar(120,100,500,170); Rectangle(120,100,500,170); Rectangle(115,95,505,175);
  Bar(120,200,500,270); Rectangle(120,200,500,270); Rectangle(115,195,505,275);
  Bar(120,300,500,370); Rectangle(120,300,500,370); Rectangle(115,295,505,375);
  SetColor(14);SetTextStyle(1,0,3);
  OutTextXY(175,120,' 1) DATI PRESTABILITI ');
  OutTextXY(175,220,' 2) INPUT MANUALI ');
  OutTextXY(175,320,' 3) DATI GENERICI ');

  repeat
   SetColor(15);OutTextXY(150,410,' Scegliere un tasto da 1 a 3 '); scelta:=readkey;
  until (scelta='1') or (scelta='2') or (scelta='3');

  CloseGraph;
  RestoreCrtMode; ColoriSchermo;

  case scelta of
    '1': Automatico;
    '2': Avvio;
    '3': Generico;
  end;
End;


PROCEDURE Menu;
 var gd,gm:integer;
     scelta,ch:char;
BEGIN
  gd:=detect;initgraph(gd,gm,'');
  if graphresult <> grOk then halt(1);
  ClearViewPort;SetViewPort(0,0,GetMaxX,GetMaxY,ClipOn);
  SetFillStyle(7,7);
  Bar(120,5,500,55); Rectangle(120,5,500,55); Rectangle(115,0,505,60);
  SetCOLOR(14); SetTextStyle(1,0,5); OutTextXY(140,0,'Menu'' principale ');
  SetBkColor(9);SetColor(15);
  Setfillstyle(7,13);
  Bar(120,100,500,170); Rectangle(120,100,500,170); Rectangle(115,95,505,175);
  Bar(120,200,500,270); Rectangle(120,200,500,270); Rectangle(115,195,505,275);
  Bar(120,300,500,370); Rectangle(120,300,500,370); Rectangle(115,295,505,375);
  SetColor(14);SetTextStyle(1,0,2);
  OutTextXY(120,120,' 1) INTRODUZIONE AL PROBLEMA ');
  OutTextXY(120,220,' 2) SIMULAZIONE DEL TRASFERIMENTO ');
  OutTextXY(120,320,' 3) USCITA DAL PROGRAMMA ');

  repeat
  SetColor(15);OutTextXY(150,410,' Scegliere un tasto da 1 a 3 '); scelta:=readkey;
  until (scelta='1') OR (scelta='2') OR (scelta='3');

  CloseGraph;
  RestoreCrtMode; ColoriSchermo;

  case scelta of
    '1': Introduzione;
    '2': Simulazione;
    '3': halt;
  end;
End;

{---------------------- INIZIO Programma ------------------------------}
  Begin{main}
   Copertina;
       repeat
           Menu;
           scelta:=readkey;
       until scelta='3';
    closegraph;
End.

  {FINE Programma}
