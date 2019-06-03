Program GILERA;   {Tempo di intercettazione pianeta arrivo-sonda,
                   in un volo interplanetario alla Hohmann }

  Uses  Crt;
   Const  rad=pi/180; kgauss=0.01720209895;
          uv0=29.78469169; ua=149597870;

   Var  ft: Text;
        gg: Integer;
        a,e,p,n,r1,r2, Eps, rinflar,
        t, teta, Tp, Anomedia, MP0,MPd,MPa, rsonda,
        np,alfa0, W,Wpri,anomecc,tg2,passo,xD,yD,xsonda,ysonda,xp2,yp2 : Real;

{--------------------------- FUNZIONI --------------------------------}

      Function ArcoTANG(num,den:Real): Real;
    Var atan: Real;
   Begin
     atan:=arctan(num/den);
      if den<0 then atan:=atan+pi;
      if atan<0 then atan:=atan+2*pi;
     ArcoTANG:=atan;
  End;

      Function ArcoCOS(x:Real): Real;
   Begin
     ArcoCOS:=pi/2-arctan(x/SQRT(1-x*x));
  End;

  Procedure ColoriSchermo;
   Begin
     textcolor(14);textbackground(1);clrscr;
  End;

   Procedure ParametriEllisse;
 Begin
   a:=(r1+r2)/2;
   e:=(r2-r1)/(r1+r2);
   p:=a*(1-e*e);
   n:=kgauss/(rad*a*sqrt(a)); {Moto medio sondaø/giorno}
   Tp:=(2*pi/kgauss)*a*sqrt(a); {Periodo di rivoluzione in giorni}
 End;

   Procedure TrePianeti;
 Begin
    np:=kgauss/(rad*r2*sqrt(r2)); {Moto medio Pianeti ø/giorno}
    alfa0:=180-np*Tp/2;           {à0= Anom.Pianeti al tempo t=0}
    rinflar:=rinflar/ua;   {raggio sf.d'influenza in UA}
    Eps:=rinflar/(r2*rad);  {anomalia Pianeti in gradi}
 End;

   Procedure KEPLER(m,e:Real; VAR v:Real);
   Var variaz,tgv: Real;
   Begin
   anomecc:= m;  {Valore iniziale di E assunto pari all'anom.media M}
 repeat
    variaz:= (m + e * SIN(anomecc) - anomecc) / (1 - e * COS(anomecc));
    anomecc:= anomecc + variaz;
 UNTIL ABS(variaz) < 1e-10;
   tgv:= SIN(anomecc / 2)/COS(anomecc / 2) * SQRT((1 + e) / (1 - e));
   v:= 2 * ArcoTANG(SQRT((1 + e) / (1 - e))*SIN(anomecc / 2), COS(anomecc / 2));
   End;

{---------------------- INIZIO Programma ------------------------------}
   Begin
{   ASSIGN (ft,'GILERA.TXT');
   Rewrite(ft);                   }

     ColoriSchermo;           {gg= intervallo di tempo in giorni}

   { -------------------   Dati Marte --------------------------
            r1:=1.00; r2:=1.524; rinflar:=577230; gg:=1;  }

   { -------------------   Dati Giove -------------------------- }
            r1:=1.00; r2:=5.2; rinflar:=48207568; gg:=10;

   { -------------------   Dati Saturno --------------------------
            r1:=1.00; r2:=9.56; rinflar:=54769450; gg:=10;  }


     ParametriEllisse;
     TrePianeti;

                     {Processo iterativo}
        t:=Tp/2-60/np; {Valore iniziale t0: pianeta a 60ø prima del rendez-vous}
        Passo:=gg; {intervallo giorni}
 repeat
        MPa:=alfa0+Eps+np*t;   {Anomalia Punto Pa}
        Anomedia:= t*2*pi/Tp; {anom. Media sonda in rad}
        KEPLER(Anomedia,e,teta);   {é=Anom.VERA Sonda}
         W:=teta/rad-MPa;

            if W<0 then begin
               t:=t-passo;
               passo:=passo/10;
            end;
        writeln('  t=',t:8:4,'  é=',teta/rad:8:4,'   -> W=',W:10:6);
       t:=t+passo;
 until abs(W) < 1e-5;

   rsonda:=p/(1+e*COS(teta));

   xsonda:=rsonda*COS(teta);
   ysonda:=rsonda*SIN(teta);

   xp2:=r2*COS(teta-Eps);
   yp2:=r2*SIN(teta-Eps);

   xD:=rsonda*COS(teta)-r2*COS(teta-Eps*rad);
   yD:=rsonda*SIN(teta)-r2*SIN(teta-Eps*rad);

        writeln;
        writeln('  T/2=',Tp/2:9:4,' -> Coo_Elioc.[UA]: xD=',xD:11:8,'   yD=',yD:11:8);
        writeln('  xsonda=',xsonda:11:8,'   ysonda=',ysonda:11:8);
        writeln('  xp2=',xp2:11:8,'   yp2=',yp2:11:8);

 { Close(ft);  }
End.

  {FINE Programma}
