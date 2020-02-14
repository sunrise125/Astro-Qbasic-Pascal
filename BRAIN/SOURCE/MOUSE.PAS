{=======================================================================}
{                                                                       }
{                    unit:  MOUSE                                       }
{                           =====                                       }
{                                                                       }
{  File: MOUSE.PAS                                                      }
{                                                                       }
{  This unit contains the procedures and functions to control the       }
{  Mouse systems optical mouse, the Microsoft mouse and the PS/2 mouse  }
{  The graphical procedures are to be used with the Ega interface.      }
{                                                                       }
{  BEFORE USING ANY CALL, RUN PROGRAM MSMOUSE.COM                       }
{                                                                       }
{                                                     feb 88            }
{                                                                       }
{=======================================================================}


{$D-,R-,S-,V-}
UNIT MOUSE;


INTERFACE
{***********************************************************************}

USES Dos;

PROCEDURE SetMouse (Visible: boolean);
FUNCTION  ResetTextMouse (Mouse: char; Fore,Back: word): boolean;
FUNCTION  ResetGraphMouse: boolean;
PROCEDURE FindMouse (var X,Y: word; Var Left,Middle,Right: boolean);
PROCEDURE GotoMouse (X,Y: word);
PROCEDURE FindButton (Button: word; Pressed: boolean; var X,Y,Count: word);
PROCEDURE SetMouseRange (X1,Y1,X2,Y2: word);
PROCEDURE SetMouseSpeed (X,Y: word);





IMPLEMENTATION
{***********************************************************************}


CONST
  MOUSEINT = $33;                { MSMOUSE.COM driver interrupt    }

VAR
  MouseMode : integer;           { 0: tet, 1: ega high resolution }



FUNCTION ResetMouse : boolean;
{-----------------------------------------------------------------------}
{ in:  -                                                                }
{ out: func.   - TRUE: mousedriver installed                            }
{ rem: reset the mouse and test INTR 51                                 }
{                                                                       }
{      THIS FUNCTION MUST BE CALLED AFTER A CHANGE FROM TEXT TO GRAPHIX }
{      OR VICE VERSA TOO                                                }
{-----------------------------------------------------------------------}

VAR
  Regs : registers;

BEGIN
  with Regs do
    begin
      AX := 0;                       { function 0 = reset   }
      Intr (MOUSEINT,Regs);
      ResetMouse := AX = $FFFF;      { -1 = mouse installed }
    end;
END; { ResetMouse }


PROCEDURE MakeTextMouse (Ascii: char; Fore,Back: word);
{-----------------------------------------------------------------------}
{ in:  Ascii - the character to use as Mouse (1..255)                   }
{      Fore  - foreground color 0..15                                   }
{      Back  - background color 0..7 or 8..15 for blinking              }
{ out: -                                                                }
{ rem: define the mouse pattern to display                              }
{-----------------------------------------------------------------------}

VAR
  Regs       : registers;

BEGIN
  with Regs do
    begin
      AX := 10;   { function 10 = set text cursor }
      BX := 0;    { select software picture       }
      CX := 0;    { screen mask                   }
      DH := (Back SHL 4) OR Fore;    { attribute  }
      DL := ord (Ascii);
      Intr (MOUSEINT,Regs);
    end;
END; { MakeTextMouse }



PROCEDURE ConvertMouse (PixelToChar: boolean; var X,Y: word);
{-----------------------------------------------------------------------}
{ in:  PixelToChar - TRUE: convert pixels to characters                 }
{      X,Y         - pixel (or char) column and row                     }
{ out: X,Y         - char (or pixel) column and row                     }
{ rem: convert screen coordinates to row and column or vice versa       }
{-----------------------------------------------------------------------}

BEGIN
  if MouseMode = 0 then
    if PixelToChar then
      begin
        X := X DIV 8 +1;
        Y := Y DIV 8 +1;
      end
    else
      begin
        X := (X -1) * 8;
        Y := (Y -1) * 8;
      end;
END; { ConvertMouse }


PROCEDURE SetMouse (Visible: boolean);
{-----------------------------------------------------------------------}
{ in:  Visible - TRUE: show mouse, FALSE: hide mouse                    }
{ out:                                                                  }
{ rem: show or hide the mouse                                           }
{-----------------------------------------------------------------------}

VAR
  Regs : registers;

BEGIN
  with Regs do
    begin
      if Visible then
        AX := 1               { function 1 = show mouse }
      else
        AX := 2;              { function 2 = hide mouse }
      Intr (MOUSEINT,Regs);
    end;
END; { SetMouse }



FUNCTION ResetTextMouse (Mouse: char; Fore,Back: word): boolean;
{-----------------------------------------------------------------------}
{ in:  Mouse - any character (1..255) to be used as the mouse           }
{      Fore  - foreground color of the mouse                            }
{      Back  - background color of the mouse                            }
{ out: func. - TRUE: mousedriver installed                              }
{ rem: reset mouse, make a mouse symbol                                 }
{-----------------------------------------------------------------------}

VAR
  Ok : boolean;

BEGIN
  Ok := ResetMouse;
  if Ok then
    begin
      MakeTextMouse (Mouse,Fore,Back);
      MouseMode := 0;
      ResetTextMouse := Ok;
    end
  else
    ResetTextMouse := FALSE;
END; { ResetTextMouse }




FUNCTION ResetGraphMouse: boolean;
{-----------------------------------------------------------------------}
{ in:  -                                                                }
{ out: func. - TRUE: mousedriver installed                              }
{ rem: reset mouse, make white uparrow visible                          }
{-----------------------------------------------------------------------}

VAR
  Ok : boolean;

BEGIN
  Ok := ResetMouse;
  MouseMode := 1;
  ResetGraphMouse := Ok;
END; { ResetGraphMouse }



PROCEDURE FindMouse (var X,Y: word; Var Left,Middle,Right: boolean);
{-----------------------------------------------------------------------}
{ in:  -                                                                }
{ out: X      - X position of mouse (1..80  in text       mode)         }
{                                   (0..639 if in grafics mode)         }
{      Y      - Y position of mouse (1..25  in text       mode)         }
{                                   (0..349 if in grafics mode)         }
{      Left   - TRUE if left button is down, FALSE if button is up      }
{      Middle -         middle                                          }
{      Right  -         right                                           }
{ rem: get the mouse position and the status of the push butons         }
{-----------------------------------------------------------------------}

VAR
  Regs : registers;

BEGIN
  with Regs do
    begin
      AX := 3;              { function 3: get mouse and buttons }
      Intr (MOUSEINT,Regs);
      X := CX;
      Y := DX;
      Left   := BX AND 1 = 1;
      Right  := BX AND 2 = 2;
      Middle := BX AND 4 = 4;
    end;
  ConvertMouse (TRUE,X,Y);
END; { FindMouse }


PROCEDURE GotoMouse (X,Y: word);
{-----------------------------------------------------------------------}
{ in:  X      - X position of mouse (1..80  in text       mode)         }
{                                   (0..639 if in grafics mode)         }
{      Y      - Y position of mouse (1..25  in text       mode)         }
{                                   (0..349 if in grafics mode)         }
{ out: -                                                                }
{ rem: set the mouse position                                           }
{-----------------------------------------------------------------------}

VAR
  Regs : registers;

BEGIN
  ConvertMouse (FALSE,X,Y);
  with Regs do
    begin
      AX := 4;              { function 4 : set mouse }
      CX := X;
      DX := Y;
      Intr (MOUSEINT,Regs);
    end;
END; { GotoMouse }



PROCEDURE FindButton (Button: word; Pressed: boolean; var X,Y,Count: word);
{-----------------------------------------------------------------------}
{ in:  Button  - 0: left   1: right   2: middle                         }
{      Pressed - TRUE: button press,  FALSE: button release             }
{ out: X       - X position of mouse (1..80  in text       mode)        }
{                                    (0..639 if in grafics mode)        }
{      Y       - Y position of mouse (1..25  in text       mode)        }
{                                    (0..349 if in grafics mode)        }
{      Count   - count of button presses since last call                }
{ rem: get the position where the button is pressed the last time, or   }
{      released the last time                                           }
{-----------------------------------------------------------------------}

VAR
  Regs : registers;

BEGIN
  with Regs do
    begin
      if Pressed then
        AX := 5             { function 5: get button pressed  }
      else
        AX := 6;            { function 6: get button released }
      BX := Button;
      Intr (MOUSEINT,Regs);
      Count := BX;
      X := CX;
      Y := DX;
    end;
  ConvertMouse (TRUE,X,Y);
END; { FindButton }



PROCEDURE SetMouseRange (X1,Y1,X2,Y2: word);
{-----------------------------------------------------------------------}
{ in:  X1,Y1 - upper left  coordinates (0..639,0..349) or (1..80,1..25) }
{      X2,Y2 - lower right coordinates (0..639,0..349) or (1..80,1..25) }
{ out: -                                                                }
{ rem: describe the box in which the mouse has to stay                  }
{      If an EGA/VGA is used in 43/50 line mode, this routine must be   }
{      used with SetMouseRange (1,1,80,43)  or 50)                      }
{                       0..349 if in high resolution grafics mode       }
{-----------------------------------------------------------------------}

VAR
  Regs : registers;

BEGIN
  ConvertMouse (FALSE,X1,Y1);
  ConvertMouse (FALSE,X2,Y2);
  with Regs do
    begin
      AX := 7;                   { function 7: set mouse X-range }
      CX := X1;
      DX := X2;
      Intr (MOUSEINT,Regs);
      AX := 8;                   { function 8: set mouse Y-range }
      CX := Y1;
      DX := Y2;
      Intr (MOUSEINT,Regs);
    end;
END; { SetMouseRange }


PROCEDURE SetMouseSpeed (X,Y: word);
{-----------------------------------------------------------------------}
{ in:  X - sensitivity in X direction, value 1..MAXINT or 0 = default   }
{      Y - sensitivity in Y direction, value 1..MAXINT or 0 = default   }
{ out: -                                                                }
{ rem: set the sensitivity. The larger the number, the greater the      }
{      distance the mouse has to be moved. defaultvalues are 8          }
{-----------------------------------------------------------------------}

VAR
  Regs : registers;

BEGIN
  if X = 0 then X := 8;
  if Y = 0 then Y := 8;
  with Regs do
    begin
      AX := 15;            { function 15: set mouse to pixel ratio }
      CX := X;
      DX := Y;
      Intr (MOUSEINT,Regs);
    end;
END; { SetMouseSpeed }


END.
