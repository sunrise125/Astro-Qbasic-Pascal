Program Best_Sats;     { 29.05.2019 SATENG.PAS, testo in inglese

                         15.02.2006 --> Tutti i satelliti terrestri
                                        ISTANTE x ISTANTE
                         Sono utilizzati i file Dati del sito di Celestrak
                         NOME.txt. N.ro max Sats settato a 352 (lim_max memo
                         progr. Pascal) in quanto i geostazionari attuali
                         sono 344. Gli altri tipi
                         di satelliti sono al massimo una trentina, come
                         gli importanti GPS (gps-ops.txt)
     Nota IMPORTANTE 
    Il file dati delle coordinate
    della Stazione osservativa
    (coord.dat) DEVE stare nella
    stessa cartella dell'eseguibile

                       }
{$N+}
  Uses CRT,SGP_Intf,
       SGP_Init,SGP_Conv,
       SGP_Math,SGP_Time,
       SGP4SDP4,SGP_Obs;

  Const rad=pi/180;
var
  ft,fl                : Text;
  line, w, band        : String;
  g,aux4               : Real;
  satnumber,i,j,k,TotSats,mm,yy   : Integer;
  time,tsince,jdrif : double;
  pos,vel,geodetic,obs_set : vector;


  Procedure ColoriSchermo;
   Begin
     textcolor(14);textbackground(1);clrscr;
  End;

  Procedure Titolo;
   Begin;
     gotoxy(17,2); writeln('浜様様様様様様様様様様様様様様様様様様様様様様�');
     gotoxy(17,3); writeln('�  PERTURBED Motion of Terrestrial SATELLITES �');
     gotoxy(17,4); writeln('�       (High precision Ephemerides)          �');
     gotoxy(17,5); writeln('藩様様様様様様様様様様様様様様様様様様様様様様�');

     gotoxy(17,6); writeln('浜様様様様様様様様様様様様様様様様様様様様様様�');
     gotoxy(17,7); writeln('�                                             �');
     gotoxy(17,8); writeln('�                                             �');
     gotoxy(17,9); writeln('�                                             �');
    gotoxy(17,10); writeln('�                                             �');
    gotoxy(17,11); writeln('藩様様様様様様様様様様様様様様様様様様様様様様�');

 textcolor(13);
 gotoxy(27,8); writeln('(TLE format - 3 rows)');
 textcolor(11);
 gotoxy(21,9); writeln('URL: http://celestrak.com/NORAD/elements');
 textcolor(12);
 gotoxy(22,10); writeln('Station Coords inside the file: coord.dat');
 textcolor(12);
 gotoxy(18,7); write(' Name Group Sats (xxxxxxxx.ext): '); textcolor(14);
 read (w);
  End;

  Procedure CoordinateStazione;
  Var
   fll: Text;
   aux1,aux2,aux3: double;
   Begin
     ASSIGN(fll,'coord.dat');
          Reset(fll);
     readln(fll,aux1);
     readln(fll,aux2);
     readln(fll,aux3);
  geodetic[1] := aux1*rad;   {lat in rad}
  geodetic[2] := aux2*rad;   {lon in rad}
  geodetic[3] := aux3/1000;  {alt in km}
    Close(fll);
  End;

  Procedure ContaSats;
  Label Skip;
   Begin
     ASSIGN(fl,w);
          Reset(fl);
    while not eof(fl) do begin
       readln(fl,line);
          if line='' then goto Skip;  {Chiude lo screening se ci sono delle
                                       righe vuote in fondo al Dbase}
       inc(k);
    end; {while}
Skip:
  TotSats:=trunc(int(k/3));
    Close(fl);
  End;

  Procedure SLICE;
  Begin
     ASSIGN(fl,w);
          Reset(fl);
     for i:= 1 to TotSats do begin
       readln(fl,sat_name[i]);
       readln(fl,sat_data[i,1]);
       readln(fl,sat_data[i,2]);
     end; {for i}
    Close(fl);
  End;

  Procedure DataEffemeride;
 Var m,ha,mn: WORD;
        y: INTEGER;
  a,b,kk,c,d: REAL;
  Begin
   gotoxy(17,12); writeln('浜様様様様様様様様様様様様様様様様様様様様様様�');
   gotoxy(17,13); writeln('�           Ephemeris Date and Time           �');
   gotoxy(17,14); writeln('�                   (in UTC)                  �');
   gotoxy(17,15); writeln('�                                             �');
   gotoxy(17,16); writeln('�                                             �');
   gotoxy(17,17); writeln('�                                             �');
   gotoxy(17,18); writeln('�         (fractional days also accepted      �');
   gotoxy(17,19); writeln('藩様様様様様様様様様様様様様様様様様様様様様様�');

repeat gotoxy(22,15);write('     Day = ');readln(g);
      until (g<32) and (g>=0);
repeat gotoxy(22,16);write('   Month = ');readln(m);
      until m in [1..12];
       gotoxy(22,17);write('    Year = ');readln(y);
repeat gotoxy(43,15);write(' Hour UT = ');readln(ha);
      until ha in [0..23];
repeat gotoxy(43,16);write(' Min. UT = ');readln(mn);
      until mn in [0..59];

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
      jdrif:=b+c+d+g+1720994.5+ha/24+mn/1440;
  aux4:=g+ha/24+mn/1440;
  End;


{------------------------ INIZIO PROGRAMMA ----------------}
BEGIN
  ColoriSchermo;
  Titolo;
  CoordinateStazione;    {.... OK verificata l'assegnazione delle CooStaz.}
  ContaSats;
  SLICE;
  DataEffemeride;

  ClrScr;
    writeln(' Observer Coords. (lat,lon,H)= ', geodetic[1]/rad:6:4,'  ',geodetic[2]/rad:6:4,'  ',geodetic[3]*1000:4:1);
    for satnumber := 1 to TotSats do
    begin
    Writeln(sat_name[satnumber]);
    Writeln(sat_data[satnumber,1]);
    Writeln(sat_data[satnumber,2]);

    Convert_Satellite_Data(satnumber);
 time := Julian_Date_of_Epoch(epoch);
 tsince := (jdrif-time)*1440;   {t-span in Minuti da Epoca TLE}
 {------------------------------------------------}
      if ideep = 0 then
        begin
        {Write('SGP4');}
        SGP4(tsince,iflag,pos,vel);
        end {if SGP4}
      else
        begin
        {Write('SDP4');}
        SDP4(tsince,iflag,pos,vel);
        end; {else SDP4}
 {------------------------------------------------}
      Convert_Sat_State(pos,vel);
      Calculate_Obs(pos,vel,geodetic,jdrif,obs_set);

      writeln(' JD0= ',time:15:8,'    dt= ',tsince/1440:10:8,'   ideep= ',ideep);
      textcolor(10);
      write(' ',yy:4,'/',mm:1,'/',aux4:8:5);
      textcolor(14);
      write(' Azim.[deg]= ',obs_set[1]/rad:8:4,'   Height[deg]= ',obs_set[2]/rad:8:4);
          if obs_set[2] < 0 then begin
               textcolor(12);
               band:='   NOT-Visible';
               end
             else begin
               textcolor(11);
             band:='   Visible';
          end;
      writeln(band); textcolor(14);
      writeln(' Topoc.Dist.(`range`) [km]= ',obs_set[3]:9:3);
      Write(' (X,Y,Z)={',pos[1]:8:0,pos[2]:8:0,pos[3]:8:0,'} -> ');
      Writeln('r_geoc= ',sqrt(pos[1]*pos[1]+pos[2]*pos[2]+pos[3]*pos[3]):6:0,' km');
      writeln('------------------');
      {repeat until keypressed; readln;}
    end; {for satnumber}

END.
{-------------- FINE PROGRAMMA ------------}
