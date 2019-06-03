 Program DEC_TYC;    {DECoding TYCho Catalogues}

  Uses  Crt,Dos;
   Var
      fb          : file of Byte;
      fascii      : Text;
      k,j,SUP,INF,nro,mv,RA,DE   : LongInt;
      byt         : Byte;
      nro1,nro2,nro3,nro4, ar1,ar2,ar3, de1,de2,de3  :Word;
      intero,frazio,RAgra,DEgrad : Real;
      nr_regio,nr_star, nr_extra, number: Integer;
      nr32: Longint;
      SegnoDE: Char;
      {---------------------------------------------------------------}
      ARint,ARore,ARmin,DEgra,DEpri :Integer;
      ARsec,ARminfraz,DEsec,DEprifraz: Real;

 { STARTING Program}
   Begin
    SegnoDE:=' ';
     textcolor(14);textbackground(1);clrscr;
   gotoxy(22,1); writeln('Binary file reading STARTED');

   ASSIGN (fb,'c:\Temp\hnsky2\tyc_h002.d32');    {BINARY original file to be decoded}
   Reset(fb);

   ASSIGN (fascii,'TYC02.TXT');                {OUTPUT final file}
   Rewrite(fascii);

   INF:=111;    {First record}
   SUP:=144417;    { Last record  .... [filesize(fb)-111]/11 - 1
                           calcolo manuale filesize si legge con dir nomefile }

   For k:=INF to SUP do begin
    For j:=1 to 11 do begin
     seek(fb,111+11*(k-1)+j-1);       {111+... because binary Dbase
                                       starts at byte 112}
     read(fb,byt);
          if j=1 then nro1:=ord(byt);   {Cat.nr}
          if j=2 then nro2:=ord(byt);
          if j=3 then nro3:=ord(byt);
          if j=4 then nro4:=ord(byt);

          if j=5 then ar1:=ord(byt);    {Right Ascension}
          if j=6 then ar2:=ord(byt);
          if j=7 then ar3:=ord(byt);

          if j=8  then de1:=ord(byt);   {Declination}
          if j=9  then de2:=ord(byt);
          if j=10 then de3:=ord(byt);

          if j=11 then mv:=ord(byt);    {Magnitude}
    end; {...for j= ...}

    nr32:=256*256*256*nro4+256*256*nro3+256*nro2+nro1;

    nr_extra:=1;
    nr_regio:=((-nr32) and $3FFF0000) shr 16;
    nr_star:=(-nr32) and $7FFF;
    if (((-nr32) and $40008000)>0) then {extra tycho number}
    begin
    if (((-nr32) and $40000000)>0) then nr_extra:=2
      else nr_extra:=3; {HEX-00008000 bit set}
    end;

      RA:=256*256*ar3+256*ar2+ar1;
      RAgra:=RA/(256*256*256-1)*360;

        ARore:=trunc(RAgra/15);
        ARminfraz:=frac(RAgra/15)*60;
        ARmin:=trunc(ARminfraz);
        ARsec:=frac(ARminfraz)*60;

  {---------------------------------------------------------------}
      DE:=256*256*de3+256*de2+de1;
           if de3>=128 then DE:=DE-256*256*256;
      DEgrad:=DE/(256*256*256-2)*180;
           if de3<=128 then SegnoDE:='+' else SegnoDE:='-';

        DEgra:=trunc(abs(DEgrad));
        DEprifraz:=frac(abs(DEgrad))*60;
        DEpri:=trunc(DEprifraz);
        DEsec:=frac(DEprifraz)*60;
  {---------------------------------------------------------------}

      if mv>=128 then mv:=mv-256;
      if mv<0 then mv:=mv+256;

 write(k:6,')  ',nr_regio:1,'  ',nr_star:4,' ',nr_extra:1,'  ',ARore:4,ARmin:3,ARsec:7:3);
 writeln('   ',+SegnoDE,DEgra:2,DEpri:3,DEsec:6:2,'   ',mv/10:4:1);

 write(fascii,k:6,')  ',nr_regio:1,'  ',nr_star:4,' ',nr_extra:1,'  ',ARore:4,ARmin:3,ARsec:7:3);
 writeln(fascii,'   ',+SegnoDE,DEgra:2,DEpri:3,DEsec:6:2,'   ',mv/10:4:1);

{ AR e DE in gradi
 write(k:6,')  ',nr_regio:1,'  ',nr_star:4,' ',nr_extra:1,'  ',RAgra:12:8);
 writeln('   ',DEgrad:12:8,'   ',mv/10:4:1);
}
 end; {...for k= ...}

     Close(fb);
     Close(fascii);

 writeln('                     Binary file reading ENDED');

  End.

   {ENDING Program }

