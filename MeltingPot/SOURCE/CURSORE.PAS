UNIT Cursore;
INTERFACE
{$X+}

USES Graph,Dos,Crt,Drivers,Objects,Mouse;

CONST
PathBGI='c:\tp\bgi';
MargineSx=30;         (* margini *)
MargineSopra=40;
Marg=3;
IIICol=300;
UltimaRiga=12;
ColoreBase=lightGray;
ColoreSfondo=white;    (* colori *)
ColoreDescr=Black;
BarraInput=Black;
ColRetMossa=white;
ColorePulsante=LightGray;
ColoreOmbra=DarkGray;
ColoreCornice=Blue;
ColoreRettangolo=Black;
ColStrPulsante=Black;
ColoreParola=Black;
ColoreParolaScroll=White;
FinCentrato=$0001;
FinAltoSx=$0002;
FinAltoDx=$0003;
MaxLista=100;          (* max Array *)
MaxDesc=20;
MaxCol=3;
Step=2;
Spazio=2;
SDi=20;                (*lunghezza descrizioni*)
SDe=20;
SHDe=70;
SDp=2;
SDData=3;
LDe=15;
LDP=2;
LDData=10;




TYPE
SDitta=String[SDi];
SDesc=String[SDe];
SHelpDesc=String[SHDe];
TLista=Array[1..MaxLista] of SDitta;
ArrDesc=Array[1..MaxDesc] of SDesc;
ArrHelpDesc=Array[1..MaxDesc] of SHelpDesc;
TRecAtt=Array[1..MaxCol,1..MaxDesc] of String[LDe];
TImmagine=array[0..0] of byte;
PImmagine=^TImmagine;

PDataPuls=^TDataPuls;
TDataPuls=object       { Valori per lista Pulsanti }
       D:PString;      {Descrizione}
       C:Word;         {Comando}
          end;

Ppulsante=^Tpulsante;
TPulsante=object
          R:TRect;              {dimensioni rettangolo pulsante}
          PulSCom:PString;      {Stringa Descrizione Comando}
          PulTCom:Word;         {Valore Comando kbXXX}
          PulTAlt:Word;      {Valore Alt+Tasto}
          CONSTRUCTOR Init(AR:TRect;ASCom:string;ACom,AComAlt:Word);
          DESTRUCTOR Done;VIRTUAL;
          end;

PPulsColl=^TPulsColl;
TPulsColl=object(TCollection)
 PROCEDURE FreeItem(P:Pointer);VIRTUAL;
 PROCEDURE Pulsante(P:PPulsante);
 PROCEDURE PulsantePremuto(P:PPulsante);
          end;

PPunto=^TPunto;
TPunto=Object
 X,Y:Integer;
 Alt,Lar:Integer;
 Size:Word;
 MinX,MinY:Integer;
 MaxX,MaxY:Integer;       {Max Coordinate Schermo}
 MargineSx,MargineSopra:Integer;
 MargineDxF,MargineSxF,MargineSopraF,MargineSottoF:Integer;
 MargineDxOmbra,MargineSxOmbra,MargineSopraOmbra,MargineSottoOmbra:Integer;
 ImmagineSfondo,ImmagineOmbra:PImmagine;
 OptFinestra:Word;
 ColPrec:Word;
 CONSTRUCTOR Init;
 DESTRUCTOR Done;VIRTUAL;

       End;



PCursore=^TCursore;
TCursore=object(TPunto)
 StepX,StepY:Integer;     {Salti Assi X,Y}
 MinDoveX,MinDoveY,MaxDoveX,MaxDoveY:Integer;
 DoveX,DoveY:Integer;
 FineCursore:Boolean;
 CONSTRUCTOR Init;
 DESTRUCTOR Done;Virtual;
 PROCEDURE Seleziona;
 PROCEDURE Eventi(Var E:TEvent; PulsColl:PPulsColl);VIRTUAL;
            end;

PInputStr=^TInputStr;
TInputStr=Object(TCursore)
 CONSTRUCTOR Init(MarDx,MarSotto:Integer;Colore:Word);
 DESTRUCTOR Done;Virtual;
 FUNCTION  Cursore(Stringa:String;XXs,YYs,MaxLunStr:Integer;
                   PulsColl:PPulsColl):String;
           end;

PInputNum=^TInputNum;
TInputNum=Object(TCursore)
 CONSTRUCTOR Init(MarDx,MarSotto:Integer;Colore:Word);
 DESTRUCTOR Done;Virtual;
 FUNCTION  Cursore(Stringa:String;XXs,YYs,MaxLunStr:Integer;
                   PulsColl:PPulsColl):String;
           end;

PScrollAtt=^TScrollAtt;
TScrollAtt=Object(TCursore)
 DS:ArrDesc;            (* descrizioni Scheda *)
 HDS:ArrHelpDesc;       (* help descrizioni Scheda *)
 Dati:TRecAtt;
 DFile:File of TRecAtt;
 PulsScrollAtt:PPulsColl;
 PSalva,PEsci,PCambioDitta:PPulsante;
 IColonna,IIColonna,IIIColonna:Integer;  (* Posizioni X nelle Colonne *)
 PImmagineAtt:PImmagine;
 SizeAtt:Word;
 CONSTRUCTOR Init;
 DESTRUCTOR Done;VIRTUAL;
 PROCEDURE CursoreDestra(var XX,DoveXX:Integer;StepXX:Integer);VIRTUAL;
 PROCEDURE CursoreSinistra(var XX,DoveXX:Integer;StepXX:Integer);VIRTUAL;
 PROCEDURE CursoreSopra(Var YY,DoveYY:Integer;StepYY:Integer);VIRTUAL;
 PROCEDURE CursoreSotto(Var YY,DoveYY:Integer;StepYY:Integer);VIRTUAL;
 PROCEDURE Cursore;
 PROCEDURE ScrollBar(XX,YY,Step,Colore:Integer);VIRTUAL;
 PROCEDURE Salva(NF:String;D:TRecAtt);
 PROCEDURE Carica(NF:String);
 PROCEDURE Eventi(Var E:TEvent; PulsColl:PPulsColl);VIRTUAL;
        End;


PDitta=^TDitta;
TDitta=Object(TCursore)
 POk,PEsci,PCancella:PPulsante;
 PulsDitta:PPulsColl;
 CONSTRUCTOR Init;
 DESTRUCTOR Done;VIRTUAL;
 PROCEDURE Esegui;
       end;


PMsgBox=^TMsgBox;
TMsgBox=Object(TCursore)
 POk,PCancella:PPulsante;
 PulsMsgBox:PPulsColl;
 SizeMsgBox:Word;
 CONSTRUCTOR Init(S1,S2:String);
 DESTRUCTOR Done;VIRTUAL;
 PROCEDURE Esegui;
       end;

PRegAtt=^TRegAtt;
TRegAtt=object(TObject)
 RecAtt:TRecAtt;
 CONSTRUCTOR Init(R:TRecAtt);
 CONSTRUCTOR Load(VAR S:TStream);
 PROCEDURE Store(var S:TStream);
        end;

PInizio=^TInizio;
TInizio=Object(TCursore)
 POk,PEsci,PCancella:PPulsante;
 PulsInizio:PPulsColl;
 CONSTRUCTOR Init;
 DESTRUCTOR Done;VIRTUAL;
 PROCEDURE Esegui;
       end;

CONST
RTRegAtt:TStreamRec=(
   ObjType: 100;
   VMTLink:ofs(TypeOf(TRegAtt)^);
   Load  : @TRegAtt.Load;
   Store : @TRegAtt.Store
   );

PROCEDURE RegisterObj;



VAR

AltPartenza,LarPartenza:Integer;   {Memorizza Altezza e Larghezza prima
                                                         di SetTextStyle}
{-- Uscita --}
Uscita:Boolean;                    {Uscita Oggetti o Fine Programma}
Cancella:Boolean;                  {Esce oggetto}
Fine:Boolean;                      {Conclude Sessione}

{-- Mouse --}
XM,YM:Word;                        {Coordinate mouse}
LB,MB,RB:Boolean;                  {Tasti Mouse Premuti}

MesOk:Boolean;                     {messggi errore avviso}
MesCanc:Boolean;

{-- Utility --}
PROCEDURE SData(S:String;X1,Y1:integer);
PROCEDURE Data(D:longint;X1,Y1:integer);
FUNCTION LOAD(Name:String):Pointer;
Procedure memoria;

IMPLEMENTATION

FUNCTION LOAD(name:string):Pointer;
VAR
Fi:file;
szImg:Word;
Pt:Pointer;
begin
   assign(fi,name);
   reset(fi,1);
     if ioresult<>0 then halt
     else
     begin
     szimg:=FileSize(fi);
     getMem(Pt,szimg);
     Blockread(fi,Pt^,szimg);
     Load:=Pt;
     FreeMem(Pt,szimg);
      close(fi);
     end
end;

PROCEDURE Data(D:longint;X1,Y1:integer);
var
s:string;
BEGIN
 setfillstyle(solidFill,black);
 bar(X1-10,Y1,X1+50,Y1+15);
 setColor(lightRed);
 setwritemode(Copyput);
 str(D,s);
 outtextxy(X1,Y1,s);
end;

PROCEDURE SData(S:String;X1,Y1:integer);
BEGIN
 setfillstyle(solidFill,black);
 bar(0,GetmaxY-15,GetMaxX,GetMaxY+10);
 setColor(white);
 setwritemode(Copyput);
 outtextxy(X1,GetMaxY-15,s);
 setfillStyle(solidFill,ColoreSfondo);
end;

procedure memoria;
begin
 data((MemAvail ),400,300);
 SetColor(bLACK);
end;

PROCEDURE RegisterObj;
BEGIN
 RegisterType(RTRegAtt);
END;

PROCEDURE CercaFile(NF:pathstr;P:String);
VAR
   Path:PathStr;
   D:DirStr;
   N:NameStr;
   E:ExtStr;
BEGIN
 Path:=FSearch(NF,P);
 path:=fexpand(path);
 fSplit(Path,D,N,E);
END;

PROCEDURE TPulsColl.FreeItem(P:Pointer);
BEGIN
 If P<>NIL Then Dispose(PPulsante(P),Done);
 P:=NIL;
END;

PROCEDURE TPulsColl.Pulsante(P:PPulsante);
CONST
 Marg=1;MargParola=3;
VAR
 Lar,Alt:Integer;          {Alt e Lar Parola}
 S:TextSettingsType;       {Setting di Testo}
 PosXs,PosYs:Integer;      {Posizioni per Scrittura Stringa}
 MaxPosYs:Integer;
 PX1,PX2,PY1,PY2:Integer;  {Posizionde Rettangolo}
 Stringa,StringaTilde:String;           {Descrizione Comando}
 S1,S2:String;
 Ch:Char;
 PosTilde,i:Integer;
 ColPrec:Word;
BEGIN
 FindMouse(XM,YM,LB,MB,RB);
 GotoMouse(0,0);
 S1:='';S2:='';StringaTilde:='';
 Lar:=LarPartenza;Alt:=AltPartenza;
 PX1:=(P^.R.A.X)*Lar;PY1:=(P^.R.A.Y)*Alt;
 PX2:=(P^.R.B.X)*Lar;PY2:=(P^.R.B.Y)*Alt;
 Stringa:=P^.PulSCom^;
 For i:=1 to Length(Stringa) do
  begin
   S1:=Copy(Stringa,i,1);
   if S1<>'~' then S2:=S2+S1
    else PosTilde:=i;
  end;
 Stringa:=S2;
 S1:='';
 for i:=1 to Length(Stringa) do
   if i=PosTilde then
    S1:=S1+'-'
   else
    S1:=S1+' ';
 StringaTilde:=S1;
 GetTextSettings(S);
 SetFillStyle(SolidFill,ColorePulsante);
 Bar(PX1,PY1,PX2,PY2);
 ColPrec:=GetColor;
 SetColor(Black);
 Rectangle(PX1,PY1,PX2,PY2);
 SetFillStyle(SolidFill,ColoreSfondo);
 Bar(PX1+Marg,PY1+Marg,PX1+Marg,PY2-Marg);
 Bar(PX1+Marg,PY1+Marg,PX2-Marg,PY1+Marg);
 SetFillStyle(SolidFill,DarkGray);
 Bar(PX1+Marg,PY2-Marg,PX2-Marg,PY2-Marg);
 Bar(PX2-Marg,PY1+Marg,PX2-Marg,PY2-Marg);
 PosXs:=PX1+abs((PX2-PX1) div 2);
 PosYs:=PY1+abs((PY2-PY1) div 2)-1; {-1 per aggiustare +1 di stringatilde}
 SetTextJustify(CenterText,BottomText);
 OutTextXY(PosXs,PosYs+(alt+1),StringaTilde); {+1 per abbassare ancora }
 SetTextJustify(CenterText,CenterText);
 SetColor(Black);
 OutTextXY(PosXs,PosYs,Stringa);
 SetTextJustify(S.Horiz,S.Vert);
 gotomouse(XM,YM);
 ShowMouse;
 SetColor(ColPrec);
END;

PROCEDURE TPulsColl.PulsantePremuto(P:PPulsante);
CONST Marg=1;
VAR
 Lar,Alt:Integer;          {Alt e Lar Parola}
 S:TextSettingsType;       {Setting di Testo}
 PosXs,PosYs:Integer;      {Posizioni per Scrittura Stringa}
 MaxPosYs:Integer;
 PX1,PX2,PY1,PY2:Integer;  {Posizionde Rettangolo}
 Stringa,StringaTilde:String;           {Descrizione Comando}
 S1,S2:String;
 Ch:Char;
 PosTilde,i:Integer;
 ColPrec:Word;
BEGIN
 FindMouse(XM,YM,LB,MB,RB);
 GotoMouse(0,0);
 HideMouse;
 S1:='';S2:='';StringaTilde:='';
 Lar:=LarPartenza;Alt:=AltPartenza;
 PX1:=(P^.R.A.X)*Lar;PY1:=(P^.R.A.Y)*Alt;
 PX2:=(P^.R.B.X)*Lar;PY2:=(P^.R.B.Y)*Alt;
 Stringa:=P^.PulSCom^;
 For i:=1 to Length(Stringa) do
  begin
   S1:=Copy(Stringa,i,1);
   if S1<>'~' then S2:=S2+S1
    else PosTilde:=i;
  end;
 Stringa:=S2;
 S1:='';
 for i:=1 to Length(Stringa) do
   if i=PosTilde then
    S1:=S1+'-'
   else
    S1:=S1+' ';
 StringaTilde:=S1;
 GetTextSettings(S);
 SetFillStyle(SolidFill,ColorePulsante);
 Bar(PX1,PY1,PX2,PY2);
 ColPrec:=GetColor;
 SetColor(Black);
 Rectangle(PX1,PY1,PX2,PY2);
 SetColor(DarkGray);
 SetFillStyle(SolidFill,DarkGray);
 Bar(PX1+Marg,PY1+Marg,PX1+Marg,PY2-Marg);
 Bar(PX1+Marg,PY1+Marg,PX2-Marg,PY1+Marg);
 SetFillStyle(SolidFill,ColoreSfondo);
 Bar(PX1+Marg,PY2-Marg,PX2-Marg,PY2-Marg);
 Bar(PX2-Marg,PY1+Marg,PX2-Marg,PY2-Marg);
 PosXs:=PX1+abs((PX2-PX1) div 2);
 PosYs:=PY1+abs((PY2-PY1) div 2);
 SetTextJustify(CenterText,BottomText);
 OutTextXY(PosXs,PosYs+(alt+3),StringaTilde);
 SetTextJustify(CenterText,CenterText);  {Formato Centrato}
 SetColor(Black{ColoreStr});
 OutTextXY(PosXs,PosYs,Stringa);
 SetTextJustify(S.Horiz,S.Vert);     { formato precedente a centrato }
 DELAY(300);   (* stop per xx MS *)
 GotoMouse(XM,YM);
 ShowMouse;
 SetColor(ColPrec);
END;

CONSTRUCTOR TPulsante.Init(AR:TRect;ASCom:string;ACom,AcomAlt:Word);
BEGIN
 R:=AR;
 PulSCom:=NewStr(ASCom);
 PulTCom:=ACom;
 PulTAlt:=AComAlt;
END;

DESTRUCTOR TPulsante.Done;
BEGIN
 DisposeStr(PulsCom);
END;


CONSTRUCTOR TRegAtt.Init(R:TRecAtt);
BEGIN
 RecAtt:=R;

END;

CONSTRUCTOR TRegAtt.Load(VAR S:TStream);
BEGIN
 S.Read(RecAtt,SizeOf(TRecAtt));
END;

PROCEDURE TRegAtt.Store(var S:TStream);
BEGIN
 S.Write(RecAtt,SizeOf(TRecAtt));
END;

CONSTRUCTOR TPunto.Init;
BEGIN
 X:=0;Y:=0;
 Alt:=TextHeight('M');                          {Altezza Lettere}
 Lar:= TextWidth('M');                          {Larghezza Lettere}
 SetWriteMode(CopyPut);
 OptFinestra:=FinCentrato;
 MinX:=0;MinY:=0;
 MaxX:=GetMaxX;MaxY:=GetMaxY;
 SetColor(ColoreParola);
 InitEvents;
END;

DESTRUCTOR TPunto.Done;
BEGIN

END;

CONSTRUCTOR TCursore.Init;
VAR I:Integer;
BEGIN
 TPunto.Init;
 MinX:=MargineSx;{+(SDe+Spazio)*lar;}
 MinY:=MargineSopra;
 X:=MinX;
 Y:=MinY;
 DoveX:=1;
 DoveY:=1;
 FineCursore:=False;
 Uscita:=False;
 Cancella:=False;
END;

DESTRUCTOR TCursore.Done;
BEGIN
END;


PROCEDURE TCursore.Seleziona;
BEGIN
 FineCursore:=True;
END;

PROCEDURE TCursore.Eventi(Var E:TEvent; PulsColl:PPulsColl);
var i:integer;
    PosMouse:TPoint;
    AltKey:Word;
    Pul:PPulsante;
BEGIN
  GetMouseEvent(E);
  if (E.What=evMouseDown) and (E.Buttons=mbLeftButton) then
  begin
   PosMouse:=E.Where;
   for i:=0 to (PulsColl^.Count-1) do
    begin
     Pul:=PulsColl^.At(i);
     if Pul^.R.Contains(PosMouse)
      then E.keyCode:=Pul^.PulTCom;
    end;
  end;   {Mouse}

  if E.What=evNothing then
   GetKeyEvent(E);
  if (E.What=evKeyDown) then
   begin
   AltKey:=E.KeyCode;
   for i:=0 to (PulsColl^.Count-1) do
    begin
     Pul:=PulsColl^.At(i);
     if Pul^.PulTAlt=AltKey
      then E.keyCode:=Pul^.PulTCom;
    end;
  end;  {KeyBoard}
END;

CONSTRUCTOR TInputStr.Init(MarDx,MarSotto:Integer;Colore:Word);
BEGIN
 TCursore.Init;
 MinX:=MargineSx+(SDe+Spazio)*lar;
 MinY:=MargineSopra;
 ColPrec:=GetColor;
 SetColor(Colore);
END;

DESTRUCTOR TInputStr.Done;
BEGIN
 SetColor(ColPrec);
END;

FUNCTION TInputStr.Cursore(Stringa:String;XXs,YYs,MaxLunStr:Integer;
                          PulsColl:PPulsColl):String;
VAR
 Ch:Char;
 ColRetMossa:Integer;
 XS,Ys:Integer;
 St,Sch:String;
 MinStr,MaxStr:Integer;
 E:TEvent;

FUNCTION Destra(Var S:String;PosX:Integer):String;
BEGIN
 if PosX<=Length(S) then
  Destra:=copy(S,PosX,1)
 else Destra:='';
END;

FUNCTION CancellaParola(Var S:String;PosX:Integer):String;
BEGIN
  Dec(PosX,2);
  if PosX>0 then
  CancellaParola:=copy(S,1,PosX)
 else CancellaParola:='';
END;

BEGIN

 St:='';
 Xs:=XXs;Ys:=YYs;
 MinStr:=Xs;
 MaxStr:=Xs+(MaxLunStr*lar);

 repeat
  Eventi(E,PulsColl);
  case E.KeyCode of
   KbRight:begin              {Destra}
          if Xs>=MinStr then
           begin
            SCh:=Copy(Destra(Stringa,DoveX),1,1);
            if length(Sch)>0 then
            Begin
              OuttextXY(Xs,Ys,SCh);
              St:=St+SCh;
              Inc(Xs,Lar);Inc(DoveX);
            end;
           end
         end;
   KbBack,
   KbLeft:begin              {Sinistra}
          if Xs>=MinStr then
           begin
            St:=CancellaParola(St,DoveX);
            if length(St)>=0 then
             begin
              Dec(DoveX);Dec(Xs,Lar);
              SetFillStyle(SolidFill,BarraInput);
              if  Length(St)=0 then
               begin St:='';DoveX:=1;Xs:=MinStr; end;
              Bar(Xs,Ys,Xs+Lar,Ys+alt);
             end;
           end
         end;
     KbESC:begin;
           Uscita:=True;
           St:='';
          end;
     KbUp,
     KbDown,
     KbEnter:begin
           FineCursore:=true;
          end;
     KbDel:Begin
            Cancella:=True;
           end;
  end; {Case e.Keycode}

  Case E.CharCode of
     'A'..'z',
     '0'..'9',
     '-','.','/':begin
                  Ch:=E.CharCode;
                  if Xs<MaxStr then
                   begin
                    Ch:=UpCase(Ch);
                    St:=St+Ch;
                    OutTextXY(XS,YS,Ch);
                    Inc(XS,lar);
                    Inc(DoveX);
                   end;
                 end;

    end;       {Case E.CharCode}
  until FineCursore or Uscita or Cancella;
 Cursore:=St;
END;

CONSTRUCTOR TInputNum.Init(MarDx,MarSotto:Integer;Colore:Word);
BEGIN
 TCursore.Init;
 MinX:=MargineSx+(SDe+Spazio)*lar;
 MinY:=MargineSopra;
 ColPrec:=GetColor;
 SetColor(Colore);
END;

DESTRUCTOR TInputNum.Done;
BEGIN
 SetColor(ColPrec);
END;

FUNCTION TInputNum.Cursore(Stringa:String;XXs,YYs,MaxLunStr:Integer;
                          PulsColl:PPulsColl):String;
VAR
 Ch:Char;
 ColRetMossa:Integer;
 XS,Ys:Integer;
 St,Sch:String;
 MinStr,MaxStr:Integer;
 E:TEvent;

FUNCTION Destra(Var S:String;PosX:Integer):String;
BEGIN
 if PosX<=Length(S) then
  Destra:=copy(S,PosX,1)
 else Destra:='';
END;

FUNCTION CancellaParola(Var S:String;PosX:Integer):String;
BEGIN
  Dec(PosX,2);
  if PosX>0 then
  CancellaParola:=copy(S,1,PosX)
 else CancellaParola:='';
END;

BEGIN

 St:='';
 Xs:=XXs;Ys:=YYs;
 MinStr:=Xs;
 MaxStr:=Xs+(MaxLunStr*lar);

 repeat
  Eventi(E,PulsColl);
  case E.KeyCode of
   KbRight:begin              {Destra}
          if Xs>=MinStr then
           begin
            SCh:=Copy(Destra(Stringa,DoveX),1,1);
            if length(Sch)>0 then
            Begin
              OuttextXY(Xs,Ys,SCh);
              St:=St+SCh;
              Inc(Xs,Lar);Inc(DoveX);
            end;
           end
         end;
   KbBack,
   KbLeft:begin              {Sinistra}
          if Xs>=MinStr then
           begin
            St:=CancellaParola(St,DoveX);
            if length(St)>=0 then
             begin
              Dec(DoveX);Dec(Xs,Lar);
              SetFillStyle(SolidFill,BarraInput);
              if  Length(St)=0 then
               begin St:='';DoveX:=1;Xs:=MinStr; end;
              Bar(Xs,Ys,Xs+Lar,Ys+alt);
             end;
           end
         end;
     KbESC:begin;
           Uscita:=True;
           St:='';
          end;
     KbUp,
     KbDown,
     KbEnter:begin
           FineCursore:=true;
          end;
     KbDel:Begin
            Cancella:=True;
           end;
  end; {Case e.Keycode}

  Case E.CharCode of
     '0'..'9'    :begin
                  Ch:=E.CharCode;
                  if Xs<MaxStr then
                   begin
                    Ch:=UpCase(Ch);
                    St:=St+Ch;
                    OutTextXY(XS,YS,Ch);
                    Inc(XS,lar);
                    Inc(DoveX);
                   end;
                 end
    end;       {Case E.CharCode}
  until FineCursore or Uscita or Cancella;
 Cursore:=St;
END;


CONSTRUCTOR TScrollAtt.Init;

VAR
 i,j:Integer;
 R:TRect;
 Pul:PPulsante;
BEGIN
 TCursore.Init;
 For i:=1 to MaxCol do
  for J:=1 to MaxDesc do Dati[i,j]:='';
 MinDoveX:=1;MinDoveY:=1;
 MaxDoveX:=3;MaxDoveY:=12;
 StepX:=0;
 StepY:=2*alt;
 MargineSxF:=2*lar;
 MargineDxF:=MargineSxF+(Spazio+Sde+Spazio+Lde+Spazio+Sdp+spazio+Ldp+Spazio+
             SDData+Spazio+lddata+spazio)*Lar;
 MargineSopraF:=4*AltPartenza;
 MargineSottoF:=MargineSopraF+34*AltPartenza;
 MargineSx:=5*lar;
 MargineSopra:=6*AltPartenza;
 MinX:=MargineSx+(SDe+Spazio)*lar;
 MinY:=MargineSopra;
 IColonna:=MargineSx+(Sde+Spazio)*Lar;
 IIColonna:=IColonna+(Lde+Spazio+SDP+Spazio)*Lar;
 IIIColonna:=IIColonna+(LDP+Spazio+SDData+Spazio)*Lar;
 MaxX:=IIIColonna;
 MaxY:=MinY+(UltimaRiga)*2*Alt;
 SizeAtt:=ImageSize(MargineSxF,MargineSopraF,MargineDxF,MargineSottoF);
 GetMem(PImmagineAtt,SizeAtt);
 if PImmagineAtt<>NIL then
  GetImage(MargineSxF,MargineSopraF,MargineDxF,MargineSottoF,PImmagineAtt^)
 else
  halt;
 SetFillStyle(SolidFill,ColoreSfondo);
 Bar(MargineSxF,MargineSopraF,MargineDxF,MargineSottoF);
 SetFillStyle(SolidFill,Blue);
 Bar(MargineSxF+1,MargineSopraF+1,MargineSxF+Marg-1,MargineSottoF-1); {Sx}
 Bar(MargineSxF+1,MargineSopraF+1,MargineDxF-1,MargineSopraF+Marg-1); {Su}
 Bar(MargineSxF+1,MargineSottoF-Marg+1,MargineDxF-Marg+1,MargineSottoF-1);{giu}
 Bar(MargineDxF-1,MargineSopraF+1,MargineDxF-Marg+1,MargineSottoF-1);{Dx}
 ColPrec:=GetColor;
 SetColor(ColoreRettangolo);
 Rectangle(MargineSxF,MargineSopraF,MargineDxF,MargineSottoF);
 Rectangle(MargineSxF+Marg,MargineSopraF+Marg,MargineDxF-Marg,
             MargineSottoF-Marg);
 SetColor(ColPrec);
 DS[1]:='di';
 DS[2]:='dal';
 DS[5]:='Pos.Tributaria';
 DS[6]:='Pos.CCIAA-Reg.Ditte';
 DS[7]:='Pos.Albo Artigiani';
 DS[8]:='Pos.Albo Costruttori';
 DS[9]:='Pos.Albo Autotrasp.';
 DS[10]:='R.E.C';
 DS[11]:='Lic. Comm.Minuto';
 DS[12]:='Lic. Comm.Ingrosso';
 DS[14]:='P.A. INPS';
 DS[15]:='P.A. INAIL';
 X:=MargineSx;Y:=MargineSopra;
 (* prima colonna *)
 outtextXY(X,Y,DS[5]); Inc(Y,Step*Alt); (* riga successiva *)
 outtextXY(X,Y,DS[6]); Inc(Y,Step*Alt);
 outtextXY(X,Y,DS[7]); Inc(Y,Step*Alt);
 outtextXY(X,Y,DS[8]); Inc(Y,Step*Alt);
 outtextXY(X,Y,DS[9]); Inc(Y,Step*Alt);
 outtextXY(X,Y,DS[10]);Inc(Y,Step*Alt);
 outtextXY(X,Y,DS[14]);Inc(Y,Step*Alt);
 outtextXY(X,Y,DS[14]);Inc(Y,Step*Alt);
 outtextXY(X,Y,DS[15]);Inc(Y,Step*Alt);
 outtextXY(X,Y,DS[15]);Inc(Y,Step*Alt);
 (* Licenze *)
 X:=MargineSx;
 outtextXY(X,Y,DS[11]);Inc(Y,Step*Alt);
 outtextXY(X,Y,DS[12]);Inc(Y,Step*Alt);
 (* Seconda Colonna *)
 X:=IColonna+(LDe+Spazio)*lar;
 Y:=MargineSopra;
 For i:= 1 to 12 do begin  outtextXY(X,Y,DS[1]);Inc(Y,Step*Alt); end;
 (* terza Colonna *)
 X:=IIColonna+(LDP+Spazio)*Lar;
 Y:=MargineSopra;
 For i:=1 to 12 do begin  outtextXY(X,Y,DS[2]);Inc(Y,Step*Alt); end;

 PulsScrollAtt:=New(PPulsColl,Init(10,3));
 if PulsScrollAtt<>Nil Then
  begin
   R.Assign(15,35,20,37);
   PSalva:=New(PPulsante,Init(R,'~Salva',KbF2,KbAltS));
   PulsScrollAtt^.Insert(PSalva);
   Pul:=PulsScrollAtt^.at(0);
   R.Assign(23,35,28,37);
   PEsci:=New(PPulsante,Init(R,'~Esci',KbEsc,KbAltE));
   PulsScrollAtt^.Insert(PEsci);
   Pul:=PulsScrollAtt^.at(1);
   {R.Assign(31,35,36,37);
   PCambioDitta:=New(PPulsante,Init(R,'~Ditta',KbF3,KbAltD));
   PulsScrollAtt^.Insert(PCambioDitta);
   Pul:=PulsScrollAtt^.at(2);}
   For i:=0 to (PulsScrollAtt^.Count-1) do
    begin
     Pul:=PulsScrollAtt^.At(i);
     PulsScrollAtt^.Pulsante(PulsScrollAtt^.At(i));
    end;
 end;

 X:=MinX; Y:=MinY;
 SetColor(White);
 SetFillStyle(SolidFill,BarraInput);
END;

DESTRUCTOR TScrollAtt.Done;
BEGIN
 FindMouse(XM,YM,LB,MB,RB);        {Trova e sposta Mouse}
 GotoMouse(0,0);
 putImage(MargineSxF,MargineSopraF,PImmagineAtt^,CopyPut); {Ripristina Fin}
 if PulsScrollAtt<>NIL then Dispose(PulsScrollAtt,Done);   {elimina Punt}
 if PImmagineAtt<>NIL then FreeMem(PImmagineAtt,SizeAtt);
 Uscita:=False;
 GotoMouse(XM,YM);                                         {Ripristina Mouse}
END;

PROCEDURE TScrollAtt.CursoreDestra(var XX,DoveXX:integer;StepXX:Integer);

BEGIN
 Inc(DoveXX);
 if  DoveXX>MaxDoveX then
  begin
   DoveXX:=MinDoveX;
   XX:=MinX;
  end
 else inc(XX,StepXX);

END;

PROCEDURE TScrollAtt.CursoreSinistra(var XX,DoveXX:integer;StepXX:Integer);

BEGIN
 Dec(DoveXX);
 if  DoveXX<MinDoveX then
  begin
   DoveXX:=MaxDoveX;
   XX:=MaxX;
  end
 else Dec(XX,StepXX);
END;

PROCEDURE TScrollAtt.CursoreSopra(Var YY,DoveYY:Integer;StepYY:Integer);

BEGIN
 Dec(DoveYY);
 if  DoveYY<MinDoveY then
  begin
   DoveYY:=MaxDoveY;
   YY:=MaxY;
   Dec(YY,StepYY);
  end
 else Dec(YY,StepYY);
END;


PROCEDURE TScrollAtt.CursoreSotto(Var YY,DoveYY:Integer;StepYY:Integer);

BEGIN
 Inc(DoveYY);
 if  DoveYY>MaxDoveY then
  begin
   DoveYY:=MinDoveY;
   YY:=MinY;
  end
 else Inc(YY,StepYY);
END;

PROCEDURE TScrollAtt.ScrollBar(XX,YY,Step,Colore:Integer);
VAR XX1,YY1:Integer;
BEGIN
 SetFillStyle(SolidFill,Colore);
 XX1:=XX+(step*lar);YY1:=YY+alt;
 Bar(XX,YY,XX1,YY1);
 if colore=Coloresfondo then SetColor(black)
  else SetColor(White);
 OutTextXY(XX,YY,Dati[DoveX,DoveY]);
END;

PROCEDURE TScrollAtt.Eventi(Var E:TEvent; PulsColl:PPulsColl);
VAR
PosMouse:TPoint;
BEGIN
 TCursore.Eventi(E,PulsColl);
{ if E.What=evNothing then
  begin
   GetMouseEvent(E);
   if (E.What=evMouseDown) and (E.Buttons=mbLeftButton) then
    begin
     PosMouse:=E.Where;
    end;
  end;}
END;

PROCEDURE TScrollAtt.Cursore;
VAR
 Ch:Char;
 ColRetMossa:Integer;
 XS,Ys:Integer;
 St:String;
 I:PInputStr;
 LD:Integer;
 J:Integer;
 D:PDitta;
 E:TEvent;
 C:PInizio;

BEGIN
 LD:=Lde;
 FineCursore:=False;
 Uscita:=False;
 ScrollBar(X,Y,LDe,BarraInput);
 repeat
    EVENTI(E,PulsScrollAtt);
    case E.KeyCode of
     KbF2:Begin    {Salva}
          PulsScrollAtt^.PulsantePremuto(PulsScrollAtt^.At(0));
          PulsScrollAtt^.Pulsante(PulsScrollAtt^.At(0));
         End;

     KbEsc:Begin   {Uscita}
          PulsScrollAtt^.PulsantePremuto(PulsScrollAtt^.At(1));
          PulsScrollAtt^.Pulsante(PulsScrollAtt^.At(1));
          Uscita:=True;
         End;
    { KbF3:Begin    Ditta
          PulsScrollAtt^.PulsantePremuto(PulsScrollAtt^.At(2));
          PulsScrollAtt^.Pulsante(PulsScrollAtt^.At(2));
          C:=New(PInizio,Init);
          if C<>Nil Then
           Begin
            C^.Esegui;
            Dispose(C,Done);
           end;
         End;}
     KbUp :begin              {Sopra}
          ScrollBar(X,Y,LD,ColoreSfondo);
          CursoreSopra(Y,DoveY,StepY);
          ScrollBar(X,Y,Ld,BarraInput);
         end;
     KbDown:begin              {Sotto}
          ScrollBar(X,Y,Ld,ColoreSfondo);
          CursoreSotto(Y,DoveY,StepY);
          ScrollBar(X,Y,Ld,BarraInput);
         end;
     KbRight:begin              {Destra}
          ScrollBar(X,Y,Ld,ColoreSfondo);
          CursoreDestra(X,DoveX,StepX);
          Case DoveX of      (* controlla colonna *)
               1:begin
                  X:=IColonna;
                  LD:=LDe;
                 end;
               2:begin
                  X:=IIColonna;
                  LD:=LDP;
                 end;
               3:Begin
                  X:=IIIColonna;
                  LD:=LDData;
                 end;
               end;   (* Case *)
          ScrollBar(X,Y,LD,BarraInput);
         end;
     KbLeft:begin              {Sinistra}
          ScrollBar(X,Y,Ld,ColoreSfondo);
          CursoreSinistra(X,DoveX,StepX);
          Case DoveX of      (* controlla colonna *)
               1:begin
                  X:=IColonna;
                  LD:=LDe;
                 end;
               2:begin
                  X:=IIColonna;
                  LD:=LDP;
                 end;
               3:Begin
                  X:=IIIColonna;
                  LD:=LDData;
                 end;
               end;   (* Case *)
          ScrollBar(X,Y,LD,BarraInput);
         end;
     KbEnter  :begin
           I:=New(PInputStr,Init(MargineSx,MargineSopra,ColoreParolaScroll));
           if I<>NIL then
            begin
             Dati[DoveX,DoveY]:=I^.Cursore(Dati[DoveX,DoveY],X,Y,LD,PulsScrollAtt);
             Dispose(I,Done);
            end
          end
    end;       {Case}
  until FineCursore or Uscita;
END;

PROCEDURE TScrollAtt.Carica(NF:String);
BEGIN
 Assign(DFile,NF);
 (* I- *)
 Reset(DFile);
 (* I+ *)
 if IoResult=0 then
  begin
   Read(DFile,Dati);
   Close(DFile);
  end;
END;

PROCEDURE TScrollAtt.Salva(NF:String;D:TRecAtt);
BEGIN
 Assign(DFile,NF);
 (* I- *)
 Rewrite(DFile);
 (* I+ *)
 if IoResult=0 then
  begin
   Write(DFile,Dati);
   Close(DFile);
  end;
END;

CONSTRUCTOR TDitta.Init;
CONST Tab=2;
VAR BordoOmbraX:Integer;
    R:TRect;
    Pul:PPulsante;
    i:Integer;
BEGIN
 TPunto.Init;
 MargineSxF:=15*lar;
 MargineDxF:=MargineSxF+34*Lar;
 MargineSopraF:=8*Alt;
 MargineSottoF:=MargineSopraF+10*Alt;
 MargineSxOmbra:=MargineSxF+tab*Lar;
 MargineDxOmbra:=MargineDxF+tab*Lar;
 MargineSopraOmbra:=MargineSopraF+tab*Alt;
 MargineSottoOmbra:=MargineSottoF+tab*lar;
 MargineSx:=MargineSxF+2*lar;
 MargineSopra:=MargineSopraF+2*Alt;
 Size:=ImageSize(MargineSxF,MargineSopraF,MargineDxOmbra,MargineSottoOmbra);
 GetMem(ImmagineSfondo,Size);
 GetMem(ImmagineOmbra,Size);
 if (ImmagineOmbra<>NIL) and (ImmagineSfondo<>Nil) then
  Begin
   GetImage(MargineSxF,MargineSopraF,MargineDxOmbra,MargineSottoOmbra,ImmagineSfondo^);
   SetFillStyle(SolidFill,ColoreOmbra);
   Bar(MargineSxF,MargineSopraF,MargineDxOmbra,MargineSottoOmbra);
   GetImage(MargineSxOmbra,MargineSopraOmbra,MargineDxOmbra,MargineSottoOmbra,ImmagineOmbra^);
   PutImage(MargineSxF,MargineSopraF,ImmagineSfondo^,NormalPut);
   PutImage(MargineSxOmbra,MargineSopraombra,ImmagineOmbra^,XorPut);
   FreeMem(ImmagineOmbra,Size);
  end;
 SetFillStyle(SolidFill,ColoreSfondo);
 Bar(MargineSxF,MargineSopraF,MargineDxF,MargineSottoF);
 SetFillStyle(SolidFill,ColoreCornice);
 Bar(MargineSxF+1,MargineSopraF+1,MargineSxF+Marg-1,MargineSottoF-1); {Sx}
 Bar(MargineSxF+1,MargineSopraF+1,MargineDxF-1,MargineSopraF+Marg-1); {Su}
 Bar(MargineSxF+1,MargineSottoF-Marg+1,MargineDxF-Marg+1,MargineSottoF-1);{Giu}
 Bar(MargineDxF-1,MargineSopraF+1,MargineDxF-Marg+1,MargineSottoF-1);
 Rectangle(MargineSxF,MargineSopraF,MargineDxF,MargineSottoF);
 Rectangle(MargineSxF+Marg,MargineSopraF+Marg,MargineDxF-Marg,
             MargineSottoF-Marg);

 PulsDitta:=New(PPulsColl,Init(10,3));
 if PulsDitta<>Nil Then
  begin
   R.Assign(14,17,19,19);
   POk:=New(PPulsante,Init(R,'~Ok',KbEnter,KbAltO));
   PulsDitta^.Insert(POk);
   Pul:=PulsDitta^.At(0);
   R.Assign(21,17,26,19);
   PEsci:=New(PPulsante,Init(R,'~Esci',KbEsc,KbAltE));
   PulsDitta^.Insert(PEsci);
   R.Assign(28,17,33,19);
   PCancella:=New(PPulsante,Init(R,'~Canc',KbDel,KbAltC));
   PulsDitta^.Insert(PCancella);
   For i:=0 to (PulsDitta^.Count-1) do
    begin
     PulsDitta^.Pulsante(PulsDitta^.At(i));
    end;
   Pul:=PulsDitta^.AT(0);
 end;
END;

DESTRUCTOR TDitta.Done;
var pul:ppulsante;
BEGIN
 PutImage(MargineSxF,MargineSopraF,ImmagineSfondo^,CopyPut);
 if ImmagineSfondo<>Nil then Freemem(ImmagineSfondo,Size);
 if PulsDitta<>Nil then Dispose(PulsDitta,Done);
 Uscita:=False;
END;

PROCEDURE TDitta.Esegui;
CONST LD=3;
VAR I:PInputStr;
    A:string;
BEGIN
 New(I,Init(MargineSx,MargineSopra,ColoreParolaScroll));
 if I<>Nil Then
  begin
   outtextxy(MargineSx,MargineSopra+alt,'Ditta');
   SetFillStyle(SolidFill,BarraInput);
   Bar(MargineSx+(7)*Lar,MargineSopra+alt,MargineSx+
   (7+Ld)*Lar,MargineSopra+2*alt);
   SetColor(White);
   repeat
   A:=I^.Cursore('',MargineSx+(7)*Lar,MargineSopra+alt,Ld,PulsDitta);
   if Length(A)=0 then I^.FineCursore:=False;
   until (Length(A)>0) or Uscita or Cancella;
   Dispose(I,Done);
  end;
 FindMouse(XM,YM,LB,MB,RB);
 gotoMouse(0,0);
 if (Length(A)>0)  and ((Uscita=false) and (cancella=false)) then
  begin
   PulsDitta^.PulsantePremuto(PulsDitta^.At(0));
   PulsDitta^.Pulsante(PulsDitta^.At(0));
  end;

 if Uscita then
  begin
   PulsDitta^.PulsantePremuto(PulsDitta^.At(1));
   PulsDitta^.Pulsante(PulsDitta^.At(1));
  end;

 if Cancella then
  begin
   PulsDitta^.PulsantePremuto(PulsDitta^.At(2));
   PulsDitta^.Pulsante(PulsDitta^.At(2));
  end;

END;


CONSTRUCTOR TMsgBox.Init(S1,S2:String);
CONST Tab=2;
VAR BordoOmbraX:Integer;
    R:TRect;
    Pul:PPulsante;
    i:Integer;
BEGIN
 TPunto.Init;
 MargineSxF:=15*lar;
 MargineDxF:=MargineSxF+34*Lar;
 MargineSopraF:=8*Alt;
 MargineSottoF:=MargineSopraF+10*Alt;
 MargineSxOmbra:=MargineSxF+tab*Lar;
 MargineDxOmbra:=MargineDxF+tab*Lar;
 MargineSopraOmbra:=MargineSopraF+tab*Alt;
 MargineSottoOmbra:=MargineSottoF+tab*lar;
 MargineSx:=MargineSxF+2*lar;
 MargineSopra:=MargineSopraF+2*Alt;
 SizeMsgBox:=ImageSize(MargineSxF,MargineSopraF,MargineDxOmbra,MargineSottoOmbra);
 GetMem(ImmagineSfondo,SizeMsgBox);
 GetMem(ImmagineOmbra,SizeMsgBox);
 if (ImmagineOmbra<>NIL) and (ImmagineSfondo<>Nil) then
  Begin
   GetImage(MargineSxF,MargineSopraF,MargineDxOmbra,MargineSottoOmbra,ImmagineSfondo^);
   SetFillStyle(SolidFill,ColoreOmbra);
   Bar(MargineSxF,MargineSopraF,MargineDxOmbra,MargineSottoOmbra);
   GetImage(MargineSxOmbra,MargineSopraOmbra,MargineDxOmbra,MargineSottoOmbra,ImmagineOmbra^);
   PutImage(MargineSxF,MargineSopraF,ImmagineSfondo^,NormalPut);
   PutImage(MargineSxOmbra,MargineSopraombra,ImmagineOmbra^,XorPut);
   FreeMem(ImmagineOmbra,SizeMsgBox);
  end;
 SetFillStyle(SolidFill,ColoreSfondo);
 Bar(MargineSxF,MargineSopraF,MargineDxF,MargineSottoF);
 SetFillStyle(SolidFill,ColoreCornice);
 Bar(MargineSxF+1,MargineSopraF+1,MargineSxF+Marg-1,MargineSottoF-1); {Sx}
 Bar(MargineSxF+1,MargineSopraF+1,MargineDxF-1,MargineSopraF+Marg-1); {Su}
 Bar(MargineSxF+1,MargineSottoF-Marg+1,MargineDxF-Marg+1,MargineSottoF-1);{Giu}
 Bar(MargineDxF-1,MargineSopraF+1,MargineDxF-Marg+1,MargineSottoF-1);
 Rectangle(MargineSxF,MargineSopraF,MargineDxF,MargineSottoF);
 Rectangle(MargineSxF+Marg,MargineSopraF+Marg,MargineDxF-Marg,
             MargineSottoF-Marg);

 PulsMsgBox:=New(PPulsColl,Init(10,3));
 if PulsMsgBox<>Nil Then
  begin
   R.Assign(14,17,19,19);
   POk:=New(PPulsante,Init(R,'~Ok',KbEnter,KbAltO));
   PulsMsgBox^.Insert(POk);
   R.Assign(28,17,33,19);
   PCancella:=New(PPulsante,Init(R,'~Canc',KbDel,KbAltC));
   PulsMsgBox^.Insert(PCancella);
   For i:=0 to (PulsMsgBox^.Count-1) do
    begin
     PulsMsgBox^.Pulsante(PulsMsgBox^.At(i));
    end;
 end;
 MesOk:=False;
 MesCanc:=False;
END;

DESTRUCTOR TMsgBox.Done;
var pul:ppulsante;
BEGIN
 if ImmagineSfondo<>Nil then Freemem(ImmagineSfondo,SizeMsgBox);
 if PulsMsgBox<>Nil then Dispose(PulsMsgBox,Done);
 Uscita:=False;
END;

PROCEDURE TMsgBox.Esegui;
CONST LD=3;
VAR E:TEvent;
    A:string;
BEGIN
 repeat
  Eventi(E,PulsMsgBox);
  Case E.KeyCode of
     KbF2:begin
           uscita:=True;
           PulsMsgBox^.PulsantePremuto(PulsMsgBox^.At(0));
           PulsMsgBox^.Pulsante(PulsMsgBox^.At(0));
           MesOk:=True;
          end;
     KbF3:Begin
           uscita:=true;
           PulsMsgBox^.PulsantePremuto(PulsMsgBox^.At(1));
           PulsMsgBox^.Pulsante(PulsMsgBox^.At(1));
           MesCanc:=true;
          end;
  end; {case}
 until Uscita;
 FindMouse(XM,YM,LB,MB,RB);
 gotoMouse(0,0);
 if (Length(A)>0)  and ((Uscita=false) and (cancella=false)) then
  begin
   PulsMsgBox^.PulsantePremuto(PulsMsgBox^.At(0));
  end;

 if Cancella then
  begin
   PulsMsgBox^.PulsantePremuto(PulsMsgBox^.At(2));
  end;

 PutImage(MargineSxF,MargineSopraF,ImmagineSfondo^,CopyPut);
END;

CONSTRUCTOR TInizio.Init;
CONST Tab=2;
VAR BordoOmbraX:Integer;
    R:TRect;
    Pul:PPulsante;
    i:Integer;
BEGIN
 TPunto.Init;
 Lar:=LarPartenza;Alt:=AltPartenza;
 MargineSxF:=15*lar;
 MargineDxF:=MargineSxF+20*Lar;
 MargineSopraF:=6*Alt;
 MargineSottoF:=MargineSopraF+20*Alt;
 MargineSxOmbra:=MargineSxF+tab*Lar;
 MargineDxOmbra:=MargineDxF+tab*Lar;
 MargineSopraOmbra:=MargineSopraF+tab*Alt;
 MargineSottoOmbra:=MargineSottoF+tab*lar;
 MargineSx:=MargineSxF+2*lar;
 MargineSopra:=MargineSopraF+2*Alt;
 Size:=ImageSize(MargineSxF,MargineSopraF,MargineDxOmbra,MargineSottoOmbra);
 GetMem(ImmagineSfondo,Size);
 GetMem(ImmagineOmbra,Size);
 if (ImmagineOmbra<>NIL) and (ImmagineSfondo<>Nil) then
  Begin
   GetImage(MargineSxF,MargineSopraF,MargineDxOmbra,MargineSottoOmbra,ImmagineSfondo^);
   SetFillStyle(SolidFill,ColoreOmbra);
   Bar(MargineSxF,MargineSopraF,MargineDxOmbra,MargineSottoOmbra);
   GetImage(MargineSxOmbra,MargineSopraOmbra,MargineDxOmbra,MargineSottoOmbra,ImmagineOmbra^);
   PutImage(MargineSxF,MargineSopraF,ImmagineSfondo^,NormalPut);
   PutImage(MargineSxOmbra,MargineSopraombra,ImmagineOmbra^,XorPut);
   FreeMem(ImmagineOmbra,Size);
  end;
 SetFillStyle(SolidFill,ColoreSfondo);
 Bar(MargineSxF,MargineSopraF,MargineDxF,MargineSottoF);
 SetFillStyle(SolidFill,ColoreCornice);
 Bar(MargineSxF+1,MargineSopraF+1,MargineSxF+Marg-1,MargineSottoF-1); {Sx}
 Bar(MargineSxF+1,MargineSopraF+1,MargineDxF-1,MargineSopraF+Marg-1); {Su}
 Bar(MargineSxF+1,MargineSottoF-Marg+1,MargineDxF-Marg+1,MargineSottoF-1);{Giu}
 Bar(MargineDxF-1,MargineSopraF+1,MargineDxF-Marg+1,MargineSottoF-1);
 Rectangle(MargineSxF,MargineSopraF,MargineDxF,MargineSottoF);
 Rectangle(MargineSxF+Marg,MargineSopraF+Marg,MargineDxF-Marg,
             MargineSottoF-Marg);
 OutTextXY(MargineSx,MargineSopra,'1. Tabella Attivia');
 OutTextXY(MargineSx,MargineSopra+2*alt,'2. Tabella Vid. Registri');
 OutTextXY(MargineSx,MargineSopra+4*Alt,'3. Tabella Agg. Contabili');
 PulsInizio:=New(PPulsColl,Init(10,3));
 if PulsInizio<>Nil Then
  begin
   R.Assign(18,21,23,23);
   POk:=New(PPulsante,Init(R,'~Ok',KbEnter,KbAltO));
   PulsInizio^.Insert(POk);
   R.Assign(27,21,32,23);
   PEsci:=New(PPulsante,Init(R,'~Esci',KbEsc,KbAltE));
   PulsInizio^.Insert(PEsci);
   For i:=0 to (PulsInizio^.Count-1) do
    begin
     PulsInizio^.Pulsante(PulsInizio^.At(i));
    end;
 end;
END;

DESTRUCTOR TInizio.Done;
var pul:ppulsante;
BEGIN
 if ImmagineSfondo<>Nil then Freemem(ImmagineSfondo,Size);
 if PulsInizio<>Nil then Dispose(PulsInizio,Done);
 Fine:=True;
END;

PROCEDURE TInizio.Esegui;
CONST LD=1;
VAR I:PInputNum;
    SScelta:string;
    XDesc,YDesc:Integer;
    Xinp,YInp:Integer;
    NScelta,Code:Integer;
    TabAtt:PScrollAtt;
    TabVid:PScrollAtt;
    TabAgg:PscrollAtt;
    ColPr:Word;
    D:PDitta;
FUNCTION DatoNum(S:String;ADelta,BDelta:Word):BOOLEAN;
VAR
 N,Code:Integer;
BEGIN
 Val(S,N,Code);
 DatoNum:=((N>=ADelta) and (N<=BDelta))
END;

BEGIN
 Fine:=False;
 New(I,Init(MargineSx,MargineSopra,ColoreParolaScroll));
 if I<>Nil Then
  begin
   XDesc:=MargineSx;
   YDesc:=MargineSopra+8*Alt;
   XInp:=XDesc+8*Lar;
   YInp:=YDesc;
   ColPr:=GetColor;
   SetColor(ColoreParola);
   outtextxy(XDesc,YDesc,'Scelta');
   SetColor(ColPr);
   SetFillStyle(SolidFill,BarraInput);
   Bar(XInp,YInp,Xinp+LD*lar,Yinp+alt);
   SetColor(White);
   repeat
   repeat
    SScelta:=I^.Cursore('',XInp,YInp,Ld,PulsInizio);
    if Length(SScelta)=0 then I^.FineCursore:=False;
    SetFillStyle(SolidFill,BarraInput);
    Bar(XInp,YInp,Xinp+LD*lar,Yinp+alt);
   until DatoNum(SScelta,1,3) or Uscita or Cancella;
    FindMouse(XM,YM,LB,MB,RB);
    gotoMouse(0,0);
    if Uscita then
     begin
      PulsInizio^.PulsantePremuto(PulsInizio^.At(1));
      PulsInizio^.Pulsante(PulsInizio^.At(1));
     end
    else
     begin
      PulsInizio^.PulsantePremuto(PulsInizio^.At(0));
      PulsInizio^.Pulsante(PulsInizio^.At(0));
     end;
    Val(SScelta,NScelta,Code);
    if Uscita=false then
    begin
     D:=New(PDitta,init);     {Chiede Codice Ditta}
     if D<>NIL then
      begin
       repeat
        D^.Esegui;
       until cancella or uscita;
       Dispose(D,Done);
      end;
     Case NScelta of
      1:begin
         TabAtt:=New(PScrollAtt,Init);
         if TabAtt<>NIL then
          begin
           TabAtt^.Cursore;
           Dispose(TabAtt,Done);
          end;
        end;
      2:begin
        end;
      3:begin
        end;
     end;  {Case}
    end
    else
     Fine:=Uscita;
   SScelta:='';
   I^.FineCursore:=False;
   DatoNum(SScelta,1,3);
   until fine;
  Dispose(I,Done);
  end;
 PutImage(MargineSxF,MargineSopraF,ImmagineSfondo^,CopyPut);
END;

END.