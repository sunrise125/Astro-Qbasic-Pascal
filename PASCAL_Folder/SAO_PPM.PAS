 Program SAO_PPM;    {DECoding SAO e PPM Catalogues}

  Uses  Crt,Dos;
   Var
      fb          : file of Byte;
      fascii      : Text;
      k,j,SUP,INF,nro,mv,RA,DE   : LongInt;
      byt         : Byte;
      nro1,nro2,nro3,nro4, ar1,ar2,ar3, de1,de2,de3, spec,s1,s2  :Word;
      RAgra,DEgrad : Real;
      wo,wp,SegnoDE: Char;
      {---------------------------------------------------------------}
      ARint,ARore,ARmin,DEgra,DEpri :Integer;
      ARsec,ARminfraz,DEsec,DEprifraz: Real;

 { STARTING Program}
   Begin
    SegnoDE:=' ';
     textcolor(14);textbackground(1);clrscr;
   gotoxy(22,1); writeln('Binary file reading STARTED');

               {BINARY original file to be decoded}
   ASSIGN (fb,'c:\hnsky\SAO_HSKY.DAT');
  { ASSIGN (fb,'c:\hnsky\PPM_HSKY.DAT'); }
   Reset(fb);
              {OUTPUT final file}
   ASSIGN (fascii,'SAO.RIS');
  { ASSIGN (fascii,'PPM.RIS'); }
   Rewrite(fascii);

   INF:=1;    {First record}
   SUP:=10;   { Last record  .... [filesize(fb)-111]/11 - 1 }

   For k:=INF to SUP do begin
    For j:=1 to 11 do begin
     seek(fb,111+11*(k-1)+j-1);       {111+... because binary Dbase
                                       starts at byte 112}
     read(fb,byt);
          if j=1 then nro1:=ord(byt);   {Cat.nr}
          if j=2 then nro2:=ord(byt);
          if j=3 then nro3:=ord(byt);

          if j=4 then ar1:=ord(byt);   {Right Ascension}
          if j=5 then ar2:=ord(byt);
          if j=6 then ar3:=ord(byt);

          if j=7  then de1:=ord(byt);   {Declination}
          if j=8  then de2:=ord(byt);
          if j=9  then de3:=ord(byt);

          if j=10 then  mv:=ord(byt);   {Magnitude}

          if j=11 then spec:=ord(byt);  {Stellar Spectrum}

    end; {...for j= ...}

      nro:=256*256*nro3+256*nro2+nro1;

  {---------------------------------------------------------------}
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
  {---------------------------------------------------------------}
      s2:=spec mod 16; s1:=spec-s2;

         case s1 of
          0: wo:='O';  16: wo:='B';  32: wo:='A';  48: wo:='F';  64: wo:='G';
         80: wo:='K';  96: wo:='M'; 112: wo:='R'; 128: wo:='N'; 144: wo:='S';
        160: wo:='C'; 176: wo:='W'; 192: wo:='P'; 208: wo:='Q'; 224: wo:='+';
        240: wo:=' '; end;

         case s2 of
          0: wp:='0'; 1: wp:='1'; 2: wp:='2'; 3: wp:='3'; 4: wp:='4';
          5: wp:='5'; 6: wp:='6'; 7: wp:='7'; 8: wp:='8'; 9: wp:='9';
         10: wp:='A';11: wp:='B';12: wp:='C';13: wp:='E';14: wp:='+';
         15: wp:=' '; end;
  {---------------------------------------------------------------}

    writeln(k:6,')  ',nro:7,' ',ARore:4,ARmin:3,ARsec:7:3,'  '+SegnoDE,DEgra:2,DEpri:3,DEsec:6:2,mv/10:6:1,'  ',wo+wp);
 writeln(fascii,k:6,')  ',nro:7,' ',ARore:4,ARmin:3,ARsec:7:3,'  ',+SegnoDE,DEgra:2,DEpri:3,DEsec:6:2,mv/10:6:1,'  ',wo+wp);



{  RA e DE in gradi
    writeln(k:6,')  ',nro:7,'  ',RAgra:12:8,'   ',DEgrad:12:8,'  ',mv/10:4:1,'  ',wo+wp);
 writeln(fascii,k:6,')  ',nro:7,'  ',RAgra:12:8,'   ',DEgrad:12:8,'  ',mv/10:4:1,'  ',wo+wp);
}
 end; {...for k= ...}

     Close(fb);
     Close(fascii);

 writeln('                     Binary file reading ENDED');

  End.

   {ENDING Program }

