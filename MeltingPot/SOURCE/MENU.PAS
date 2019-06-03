PROGRAM Menu;
(*X+*)
USES graph,Crt,Dos,Objects,Cursore;


FUNCTION VideoMode:BOOLEAN;
VAR
GrDriver,GrMode:integer;
Font,Dim:integer;
S:string;
BEGIN
 GrDriver:=Detect;
 InitGraph(GrDriver,GrMode,PathBgi);
 case GrDriver of
  Cga,MCga:begin
            Font:=SmallFont;
            Dim:=4;
            end;
  Vga     :begin
            Font:=SmallFont;
            Dim:=4;
           end;
  Ega,Ega64,
  EgaMono  :begin
             Font:=SmallFont;
             Dim:=4;
            end;
  HercMono  :begin
              Font:=SmallFont;
              Dim:=4;
             end;
 else
  begin
   Font:=SmallFont;
   Dim:=6;
  end;
 end; {case}
 AltPartenza:=TextHeight('M');
 LarPartenza:=TextWidth('M');
 SetTextStyle(Font,HorizDir,Dim);
 VideoMode:=GraphResult = GrOK;
END;


PROCEDURE VideoInit;
CONST
 AltBarra=20;
VAR
 X,Y:Integer;
 S:TextSettingsType;
BEGIN
 GetTextSettings(S);
 X:=GetMaxX div 2;Y:=((AltBarra-Marg) div 2) ;
 SetColor(ColoreDescr);        {Definisce Colore Parole}
 SetFillStyle(SolidFill,ColoreBase);    {Prepara Lo Sfondo Con Colore Base}
 Bar(0,0,GetMaxX,GetMaxY);
 SetFillStyle(SolidFill,Blue);  {Disegna contorno Barra}
 Bar(0,0,GetMaxX,AltBarra);
 SetFillStyle(SolidFill,ColoreSfondo);   {Disegna Barra in Alto}
 Bar(Marg,Marg,GetMaxX-Marg,AltBarra-Marg);
 SetColor(Black);
 SetTextJustify(CenterText,CenterText);
 SetTextStyle(S.Font,HorizDir,5);
 OutTextXY(X,Y,'PROGRAMMA UTILITA');
 SetTextStyle(S.Font,S.Direction,S.CharSize);
 SetTextJustify(S.Horiz,S.Vert);
END;

PROCEDURE Presentazione;
CONST
PX1=230;
PY1=100;
PX=pX1+50;
PY=pY1+30;
PX2=pX1+150;
PY2=pY1+150;
bordo=3;
VAR
p:Pointer;
T:TextSettingsType;
Pxs,Pys:integer;
BEGIN
 GetTextSettings(T);
 SetFillStyle(SolidFill,ColoreBase);    {Prepara Lo Sfondo Con Colore Base}
 Bar(0,0,GetMaxX,GetMaxY);
 P:=Load('Pal.blk');
 SetFillStyle(SolidFill,White);
 Bar(px1,py1,px2,py2);
 PutImage(pX,pY,P^,XorPut);
 SetFillStyle(SolidFill,Blue);
 Bar(pX1,pY1,pX2,pY1+bordo);  {sopra}
 Bar(pX1,pY2,pX2,pY2-bordo);  {Sotto}
 Bar(pX1,pY1,pX1+Bordo,pY2);  {sx}
 Bar(pX2,pY1,pX2-Bordo,pY2);  {Dx}
 SetColor(Black);
 Rectangle(pX1-1,pY1-1,pX2+1,pY2+1);
 Rectangle(pX1+Bordo+1,pY1+Bordo+1,pX2-(bordo+1),pY2-(Bordo+1));
 SetTextJustify(centerText,BottomText);
 Pxs:=pX1+(abs(pX2-pX1) div 2);
 Pys:=pY+(abs(pY2-pY) div 2);
 outtextxy(PXs,PYs,'PROGRAMMA UTILITA');
 delay(1000);
 SetTextStyle(T.font,T.Direction,T.CharSize);
 SetTextJustify(T.Horiz,T.Vert);
 SetFillStyle(SolidFill,ColoreBase);    {Prepara Lo Sfondo Con Colore Base}
 Bar(0,0,GetMaxX,GetMaxY);
END;

VAR
 PRun:PInizio;
 C:PScrollAtt;
 a:string;
 s:integer;
 D:PDitta;
 i:Integer;
 P:Pointer;
BEGIN
 if VideoMode Then
  begin
   presentazione;
   VideoInit;
   PRun:=New(PInizio,Init);
   if PRun<>NIL then
    begin
     PRun^.Esegui;
    end;
  end
 else
  writeln('scheda Grafica Non Presente');
 restorecrtmode;
END.