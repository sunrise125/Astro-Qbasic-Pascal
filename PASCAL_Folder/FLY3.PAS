Program FLY3;   {TRASFERIMENTO ELLITTICO ELIOCENTRICO INTERPLANETARIO}

  Uses  Crt,Graph;
   Const  rad=pi/180; kgauss=0.01720209895; uv0=29.78469169;
          ua=149597870; muS=1;

   Var  nro,k: Integer;
        x,y,xc,yc, x1,y1,x2,y2, stile0123,spessore01,colore : Word;
        ch: Char;
        g,Rpa,Rar,mupa,muar,rinflpa,rinflar,semmagpa,semmagar,
        massrecpa,massrecar,
        a,e,p,r1,r2,             ka, t, Anomedia,
        Vp,Va,f1,f2,r1seg,r2seg,v1seg,psi,scarto,v2seg,
        Vcirc1,deltaV1,ro1,ro2, ro,
        Vcirc2,deltaV2,
        R0,vo1,vparc1,spinta1,vo2,vparc2,spinta2,
        e1,a1,e2,a2,beta,teta,gamma,
        AnomF1,timehyp1,AnomF2,timehyp2,timeelli, PerRiv,
        AnomE1,AnomE2   : Real;

  Function pot(x,y : Real) : Real;
  Begin
   if x=0 then
    Begin
     pot :=0;
     exit;
    End;
   if x<0 then
    Begin
    (* in questo caso y deve essere intero *)
     if y<>int(y) then
      Begin
       writeln('Operazione Impossibile');
       Halt;
      End;
     x:=-x;
     pot:=Exp(y*ln(x));
     if int(y/2)<>y/2 then pot:=-pot(x,y);
    end
    else pot := Exp(y*ln(x));
  End;

      Function MODULO(x,y:Real) : Real;
   Var modu : Real;
   Begin
     modu:=(x/y-INT(x/y))*y;
     MODULO:=modu;
  End;

      Function ArcoTANG(num,den:Real): Real;
    Var atan: Real;
   Begin
     atan:=arctan(num/den);
      if den<0 then atan:=atan+pi;
      if atan<0 then atan:=atan+2*pi;
     ArcoTANG:=atan;
  End;

      Function ArcoSIN(x:Real): Real;
   Begin
     ArcoSIN:=arctan(x/SQRT(1-x*x));
  End;

      Function ArcoCOS(x:Real): Real;
   Begin
     ArcoCOS:=pi/2-arctan(x/SQRT(1-x*x));
  End;

      Function SettCOSH(x:Real): Real;
   Begin
     SettCOSH:=LN(x+sqrt(x*x-1));
  End;

      Function SINH(x:Real): Real;
   Begin
     SINH:= 0.5*(EXP(x)-1/EXP(x));
  End;

  Procedure ColoriSchermo;
   Begin
     textcolor(14);textbackground(1);clrscr;
  End;

  Procedure Pausa;
   Begin
     textcolor(14);textbackground(4);
   gotoxy(64,9); write('Premere un tasto');
     textcolor(14);textbackground(1);
   ch:=readkey;
     clrscr;
  End;


  Procedure Titolo;
   Begin;
     textcolor(11);
     gotoxy(15,2); writeln('TRASFERIMENTO ELLITTICO ELIOCENTRICO INTERPLANETARIO');
     gotoxy(15,3); writeln('----------------------------------------------------');
     textcolor(12);
     gotoxy(15,4); writeln('                                  Veronica Matarazzo');
     textcolor(14);
  End;

   Procedure DatiPianeti;
   Begin;
     textcolor(14);textbackground(5);
gotoxy(10,11); write('ÚÄÄÄÄÄÄÄÄÄÄÄÄÄ Raggi e Parametri Gravitazionali ÄÄÄÄÄÄÄÄÄÄÄÄ¿');
gotoxy(10,12); write('³                   usati nel programma                     ³');
gotoxy(10,13); write('³                                                           ³');
gotoxy(10,14); write('³                  R [km]          æ [km^3/sec]             ³');
gotoxy(10,15); write('³ 1  Mercurio       2422               21 684               ³');
gotoxy(10,16); write('³ 2  Venere         6114              324 850               ³');
gotoxy(10,17); write('³ 3  Terra          6378              398 600               ³');
gotoxy(10,18); write('³ 4  Marte          3407               42 977               ³');
gotoxy(10,19); write('³ 5  Giove         71400          126 684 553               ³');
gotoxy(10,20); write('³ 6  Saturno       57540           37 911 978               ³');
gotoxy(10,21); write('³ 7  Urano         24071            5 786 122               ³');
gotoxy(10,22); write('³ 8  Nettuno       24750            6 832 062               ³');
gotoxy(10,23); write('ÀÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÙ');
     textcolor(14);textbackground(1);

  End;

   Procedure InputDati;
   Begin

gotoxy(10,6); write('Scelta del pianeta di partenza [1..8] -> ');readln(nro);
{nro:=3;}
   case nro of
   1: begin Rpa:= 2422;  mupa:=21684;     semmagpa:=0.387098;  massrecpa:=6023600 ; end;
   2: begin Rpa:= 6114;  mupa:=324850;    semmagpa:=0.723332;  massrecpa:=408523.7; end;
   3: begin Rpa:= 6378;  mupa:=398600;    semmagpa:=1;         massrecpa:=332946.0; end;
   4: begin Rpa:= 3407;  mupa:=42977;     semmagpa:=1.523688;  massrecpa:=3098708 ; end;
   5: begin Rpa:= 71400; mupa:=126684553; semmagpa:=5.202561;  massrecpa:=1047.3  ; end;
   6: begin Rpa:= 57540; mupa:=37911978;  semmagpa:=9.554747;  massrecpa:=3497.9  ; end;
   7: begin Rpa:= 24071; mupa:=5786122;   semmagpa:=19.218140; massrecpa:=22903   ; end;
   8: begin Rpa:= 24750; mupa:=6832062;   semmagpa:=30.109570; massrecpa:=19412   ; end;
   end;

gotoxy(10,7); write('Scelta del pianeta di arrivo   [1..8] -> ');readln(nro);
{nro:=4;}
              case nro of
  1: begin Rar:= 2422;  muar:=21684;     semmagar:=0.387098;   massrecar:=6023600 ; end;
  2: begin Rar:= 6114;  muar:=324850;    semmagar:=0.723332;   massrecar:=408523.7; end;
  3: begin Rar:= 6378;  muar:=398600;    semmagar:=1;          massrecar:=332946.0; end;
  4: begin Rar:= 3407;  muar:=42977;     semmagar:=1.523688;   massrecar:=3098708 ; end;
  5: begin Rar:= 71400; muar:=126684553; semmagar:=5.202561;   massrecar:=1047.3  ; end;
  6: begin Rar:= 57540; muar:=37911978;  semmagar:=9.554747;   massrecar:=3497.9  ; end;
  7: begin Rar:= 24071; muar:=5786122;   semmagar:=19.218140;  massrecar:=22903   ; end;
  8: begin Rar:= 24750; muar:=6832062;   semmagar:=30.109570;  massrecar:=19412   ; end;
  end;

 write('         Quota orb. di parcheggio del 1ø pianeta (in Km)= ');readln(ro1);
 write('         Quota orb. di parcheggio del 2ø pianeta (in Km)= ');read(ro2);

{ro1:=200; ro2:=400;}
  Pausa;

  End;

   Procedure SferaInfl;
 Begin
  rinflpa:=ua*semmagpa/pot(massrecpa,0.4);
  rinflar:=ua*semmagar/pot(massrecar,0.4);
 End;


   Procedure ParametriEllisse;
 Begin
   r1:= semmagpa; r2:= semmagar; { in UA }
   a:=(r1+r2)/2;
   e:=(r2-r1)/(r1+r2);
   p:=a*(1-pot(e,2));

   Vp:=sqrt(muS*((2/r1)-(1/a)))*uv0;
   Va:=sqrt(muS*((2/r2)-(1/a)))*uv0;

   f1:=(rinflpa/ua)/r1;  { anomalia del punto B in rad }
   f2:=pi-(rinflar/ua)/r2;  { anomalia del punto D in rad }

   r1seg:=p/(1+e*COS(f1));
   r2seg:=p/(1+e*COS(f2));

   v1seg:=sqrt(muS*((2/r1seg)-(1/a)))*uv0;
   psi:=arcoSIN((Vp*r1)/(r1seg*v1seg));
   scarto:=(PI/2)-psi;

   v2seg:=sqrt(muS*((2/r2seg)-(1/a)))*uv0;
 End;

   Procedure IperbolePartenza;
 Var aux: Real;
 Begin
   Vcirc1:=uv0*sqrt(muS/r1);
   deltaV1:=v1seg-Vcirc1;  {=v1infin}
 {--------------------------------------------------}
   R0:=Rpa+ro1;
   vo1:=sqrt(pot(deltaV1,2)+2*mupa/R0);
   vparc1:=sqrt(mupa/R0);
   spinta1:=vo1-vparc1;
   e1:=1+(R0*pot(deltaV1,2))/mupa;
   a1:=R0/(e1-1);
   beta:=arcoCOS(-(1/e1));
   teta:=(PI/2)+f1;
   gamma:=teta-beta;
 {--------------------------------------------------}
   aux:=(rinflpa+a1)/(a1*e1);
                 {tempo di percorrenza in secondi}
   AnomF1:=SettCOSH(aux);
   timehyp1:=sqrt(a1*a1*a1/mupa)*(e1*SINH(AnomF1)-AnomF1);
 End;

   Procedure EllisseTrasferimento;
 Var aux: Real;
  Begin
  aux:=sqrt((1-e)/(1+e));
  AnomE1:=2*ArcoTANG(aux*SIN(f1/2),COS(f1/2));
  AnomE2:=2*ArcoTANG(aux*SIN(f2/2),COS(f2/2));
       {tempo di percorrenza in giorni: 1/kgauss= unit… di tempo= 58.132 giorni}
  PerRiv:=2*pi*(1/kgauss)*sqrt(a*a*a/muS);
  timeelli:=(1/kgauss)*sqrt(a*a*a/muS)*((AnomE2-e*SIN(AnomE2))-(AnomE1-e*SIN(AnomE1)));
 End;

   Procedure IperboleArrivo;
 Var aux: Real;
  Begin
   Vcirc2:=uv0*sqrt(muS/r2);
   deltaV2:=v2seg-Vcirc2;  {=v2infin}
 {--------------------------------------------------}
   R0:=Rar+ro2;
   vo2:=sqrt(pot(deltaV2,2)+2*muar/R0);
   vparc2:=sqrt(muar/R0);
   spinta2:=vo2-vparc2;
   e2:=1+(R0*pot(deltaV2,2))/muar;
   a2:=R0/(e2-1);
   aux:=(rinflar+a2)/(a2*e2);
                 {tempo di percorrenza in secondi}
   AnomF2:=SettCOSH(aux);
   timehyp2:=sqrt(a2*a2*a2/muar)*(e2*SINH(AnomF2)-AnomF2);

 End;

  Procedure Risultati;
  Begin
     textcolor(14);textbackground(6);
     gotoxy(32,2); writeln('R I S U L T A T I');

     textcolor(14);textbackground(5);
  gotoxy(3,4); write('Iperbole di partenza');
     textcolor(14);textbackground(1);
  gotoxy(26,4);writeln('Veloc. di  uscita sonda:',vo1:7:3,' Km/s');
  gotoxy(26,5);writeln('Veloc. orb. di parcheg.:',vparc1:7:3,' Km/s');
  gotoxy(51,6);writeln('-----------');
  gotoxy(26,7);writeln('Consumo energetico    : ',spinta1:7:3,' Km/s');
  gotoxy(26,8);writeln('Velocit… Vì sfera infl: ',deltaV1:7:3,' Km/s');
  gotoxy(26,9); writeln('Inclinazione asintoto');
  gotoxy(26,10); writeln('    su asse focale (á): ',beta/rad:6:4,'ø');
  gotoxy(26,11); writeln('Angolo FASE asse foc. â ',gamma/rad:6:4,'ø');
  gotoxy(26,12); writeln('Tempo percorr.iperbole:   ',timehyp1/86400:5:3,' giorni');

     textcolor(14);textbackground(5);
  gotoxy(3,14); write('Ellisse di trasferim.');
     textcolor(14);textbackground(1);
  gotoxy(26,14); writeln('Semiasse maggiore (a) :  ',a:9:6,' UA');
  gotoxy(26,15); writeln('Eccentricit…      (e) :  ',e:9:6);
  gotoxy(26,16); writeln('Velocit… in B         :  ',V1seg:6:3,' km/s');
  gotoxy(26,17); writeln('Velocit… in D         :  ',V2seg:6:3,' km/s');
  gotoxy(26,18); writeln('Tempo percorr.ellisse : ',timeelli:7:3,' giorni');

     textcolor(14);textbackground(5);
  gotoxy(3,20); write('Iperbole di arrivo   ');
     textcolor(14);textbackground(1);
  gotoxy(26,20);writeln('Velocit… Vì sfera infl: ',deltaV2:7:3,' Km/s');
  gotoxy(26,21);writeln('Veloc.di ingresso sonda:',vo2:7:3,' Km/s');
  gotoxy(26,22);writeln('Veloc. orb. di parcheg.:',vparc2:7:3,' Km/s');
  gotoxy(26,23);writeln('Consumo energetico    : ',spinta2:7:3,' Km/s');
  gotoxy(26,24); writeln('Tempo percorr.iperbole:   ',timehyp2/86400:5:3,' giorni');
  ch:=readkey;
 End;

 Procedure Bilancio;
  Begin
     textcolor(14);textbackground(6);
  gotoxy(32,2); writeln('R I E P I L O G O');
     textcolor(14);textbackground(1);
  gotoxy(16,4);writeln('Consumo energetico Complessivo : ',spinta1+spinta2:7:3,' Km/s');
  gotoxy(16,6);writeln('Tempo totale di volo           : ',timehyp1/86400+timeelli+timehyp2/86400:5:3,' giorni');
{writeln;
writeln('  r1=',r1:8:6,'  r2=',r2:8:6,'  a=',a:8:6);
writeln;
writeln('  Scarto direzionale Vel_Terra e Vel_sonda= ',scarto:6:4,'ø');
}
  ch:=readkey;
 End;

   Procedure KEPLER(m,e:Real; VAR v:Real);
   Var anomecc,variaz,tgv: Real;
   Begin
   anomecc:= m;  {Valore iniziale di E assunto pari all'anom.media M}
 repeat
    variaz:= (m + e * SIN(anomecc) - anomecc) / (1 - e * COS(anomecc));
    anomecc:= anomecc + variaz;
 UNTIL ABS(variaz) < 1e-10;
   tgv:= SIN(anomecc / 2)/COS(anomecc / 2) * SQRT((1 + e) / (1 - e));
   v:= MODULO(2 * ArcoTANG(tgv, 1), 2 * pi);
   End;

  PROCEDURE PlotOrbita(a,e,teta:real);
     Begin
           ro:=ka*a*(1-e*e)/(1+e*cos(teta));
            x:=xc+round(ro*cos(teta));
            y:=yc-round(ro*sin(teta));
            Lineto(x,y);
    end;


   PROCEDURE RaggioVettore(a,e,teta:real;
         colore,stile0123,spessore13:word);
   Begin
     SetColor(colore);
     SetLineStyle(stile0123,4,spessore13);
     ro:=ka*a*(1-e*e)/(1+e*cos(teta*rad));
     x:=xc+round(ro*cos(teta*rad));
     y:=yc-round(ro*sin(teta*rad));
     Line(xc,yc,x,y);
              
     end;


  Procedure GRAFICOconica;
 Label vai;
 Const FattoreIngr=1.4; 
   Var
      ch,dt4   : String;
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
    SetTextStyle(1,0,3);
    xc:=GetmaxX div 2+Spost; yc:=GetmaxY div 2;

     SetFillStyle(1,12);
   {  SetTextStyle(0,0,1); Font,Direz.,Size
     OutTextXY(xc+250,yc-30,'Punto');OutTextXY(xc+232,yc-8,'Gamma'); }
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
        { Disegno Orbita TERRA}
         Circle(xc,yc,round(r1*ka));
   RaggioVettore (a,e,f1/rad,11,0,1); {a,e,anom.inø,colore,stile,spessore
                                  stile 0=cont;1=tratt  spess.norm=1,3=triplo}
         Circle(xc+round(r1*ka),yc,round(r1*f1*ka));
    {=========================================}
      SetColor(12);
        { Disegno Orbita MARTE}
         Circle(xc,yc,round(r2*ka));
         RaggioVettore (a,e,f2/rad,12,0,1);
         Circle(xc-round(r2*ka),yc,round(r2*(pi-f2)*ka));
    {=========================================}

        { Disegno Orbita SATELLITE per 0 ó t ó T/2 step 1 giorno}
      t:=0;
      SetTextStyle(4,0,5); {Font,Direz.,Size}
      SetColor(15);
      moveto(xc+round(r1*ka),yc); {1ø pixel grafico a inizio orbita}

     while t< PerRiv/2 do begin
        Anomedia:= t/PerRiv*2*pi; {anom. Media Orbita in rad}
        KEPLER(Anomedia,e,teta);
        PlotOrbita(a,e,teta);

              { SetColor(14);...tentativo di visualizzare i giorni!!
               str(t:3:0,dt4); OutTextXY(xc+150,yc+195,dt4);
               setcolor(1); TextHeight('    ');                   }

      t:=t+1; {step 1 giorno}
      end;

 {------------------------------------------------------------------}

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
           ka:=150;   {Ingrandimento INIZIALE}
         ClearDevice;
          ka:=150; Spost:=0; Spost1:=0;
          goto vai;
          end;

  until keypressed;

    closegraph;
    restoreCrtMode;
  End;
{-------------------------------------------------------------------------------
   INIZIO Programma}
   Begin
     ColoriSchermo;
     Titolo;
     DatiPianeti;
     InputDati;
      Pausa;
     SferaInfl;
     ParametriEllisse; IperbolePartenza; EllisseTrasferimento; IperboleArrivo;
      Risultati;
     ColoriSchermo;
     Bilancio;
     GraficoCONICA;
End.

   {FINE Programma}

  { Disegno Orbita SATELLITE
    PlotOrbita (a,e,15); }

      {Solo per prova: f1=20ø; f2=170ø, cioŠ à2=10ø -> sfere infl.}
   { RaggioVettore (a,e,20,11,0,1); a,e,anom.inø,colore,stile,spessore}
                                  {stile 0=cont;1=tratt  spess.norm=1,3=triplo
    RaggioVettore (a,e,170,12,0,1);}

        { Disegno Orbita TERRA
    PlotOrbita    (r1,0,11);  CIRCOLARE di raggio r1 UA}

        { Disegno Orbita MARTE
    PlotOrbita (r2,0,12);    CIRCOLARE di raggio r2 UA}
