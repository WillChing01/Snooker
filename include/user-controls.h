#ifndef USER-CONTROLS_H_INCLUDED
#define USER-CONTROLS_H_INCLUDED

const std::string _thinFontFile = "../../assets/Roboto-Thin.ttf";
const std::string _boldFontFile = "../../assets/Roboto-Bold.ttf";

const std::string _vertexShaderFile = "../../assets/vertex_shader.vert";
const std::string _fragmentShaderFile = "../../assets/fragment_shader.vert";
const std::string _userConfigFile = "../../assets/config.txt";
const std::string _userCueConfigFile = "../../assets/cueconfig.txt";

const std::string _cueFilePrefix = "../../assets/";

const std::map<std::string,sf::Keyboard::Key> default_controls=
{
{"Aim left",sf::Keyboard::Left},
{"Aim right",sf::Keyboard::Right},
{"Precise aim left",sf::Keyboard::Comma},
{"Precise aim right",sf::Keyboard::Period},
{"Increase power",sf::Keyboard::Up},
{"Decrease power",sf::Keyboard::Down},
{"Increase cue elevation",sf::Keyboard::LBracket},
{"Decrease cue elevation",sf::Keyboard::RBracket},
{"Offset cue tip up",sf::Keyboard::W},
{"Offset cue tip down",sf::Keyboard::S},
{"Offset cue tip left",sf::Keyboard::A},
{"Offset cue tip right",sf::Keyboard::D},
{"Strike cueball",sf::Keyboard::Space},
{"Move ball up",sf::Keyboard::Up},
{"Move ball down",sf::Keyboard::Down},
{"Move ball left",sf::Keyboard::Left},
{"Move ball right",sf::Keyboard::Right},
{"Place ball",sf::Keyboard::Enter},
{"Pause game",sf::Keyboard::Escape}
};

std::map<std::string,sf::Keyboard::Key> user_controls=default_controls;

#define ITEM(x) case sf::Keyboard::x : return #x;

std::string KeyToString(sf::Keyboard::Key k)
{
    switch(k)
    {
        ITEM(A);
        ITEM(B);
        ITEM(C);
        ITEM(D);
        ITEM(E);
        ITEM(F);
        ITEM(G);
        ITEM(H);
        ITEM(I);
        ITEM(J);
        ITEM(K);
        ITEM(L);
        ITEM(M);
        ITEM(N);
        ITEM(O);
        ITEM(P);
        ITEM(Q);
        ITEM(R);
        ITEM(S);
        ITEM(T);
        ITEM(U);
        ITEM(V);
        ITEM(W);
        ITEM(X);
        ITEM(Y);
        ITEM(Z);
        ITEM(Num0);
        ITEM(Num1);
        ITEM(Num2);
        ITEM(Num3);
        ITEM(Num4);
        ITEM(Num5);
        ITEM(Num6);
        ITEM(Num7);
        ITEM(Num8);
        ITEM(Num9);
        ITEM(Escape);
        ITEM(LControl);
        ITEM(LShift);
        ITEM(LAlt);
        ITEM(LSystem);
        ITEM(RControl);
        ITEM(RShift);
        ITEM(RAlt);
        ITEM(RSystem);
        ITEM(Menu);
        ITEM(LBracket);
        ITEM(RBracket);
        ITEM(Semicolon);
        ITEM(Comma);
        ITEM(Period);
        ITEM(Quote);
        ITEM(Slash);
        ITEM(Backslash);
        ITEM(Tilde);
        ITEM(Equal);
        ITEM(Hyphen);
        ITEM(Space);
        ITEM(Enter);
        ITEM(Backspace);
        ITEM(Tab);
        ITEM(PageUp);
        ITEM(PageDown);
        ITEM(End);
        ITEM(Home);
        ITEM(Insert);
        ITEM(Delete);
        ITEM(Add);
        ITEM(Subtract);
        ITEM(Multiply);
        ITEM(Divide);
        ITEM(Left);
        ITEM(Right);
        ITEM(Up);
        ITEM(Down);
        ITEM(Numpad0);
        ITEM(Numpad1);
        ITEM(Numpad2);
        ITEM(Numpad3);
        ITEM(Numpad4);
        ITEM(Numpad5);
        ITEM(Numpad6);
        ITEM(Numpad7);
        ITEM(Numpad8);
        ITEM(Numpad9);
        ITEM(F1);
        ITEM(F2);
        ITEM(F3);
        ITEM(F4);
        ITEM(F5);
        ITEM(F6);
        ITEM(F7);
        ITEM(F8);
        ITEM(F9);
        ITEM(F10);
        ITEM(F11);
        ITEM(F12);
        ITEM(F13);
        ITEM(F14);
        ITEM(F15);
        ITEM(Pause);
        ITEM(KeyCount);
        ITEM(Unknown);
    }
}

sf::Keyboard::Key StringToKey(std::string s)
{
    if (s=="A") {return sf::Keyboard::A;}
    else if (s=="B") {return sf::Keyboard::B;}
    else if (s=="C") {return sf::Keyboard::C;}
    else if (s=="D") {return sf::Keyboard::D;}
    else if (s=="E") {return sf::Keyboard::E;}
    else if (s=="F") {return sf::Keyboard::F;}
    else if (s=="G") {return sf::Keyboard::G;}
    else if (s=="H") {return sf::Keyboard::H;}
    else if (s=="I") {return sf::Keyboard::I;}
    else if (s=="J") {return sf::Keyboard::J;}
    else if (s=="K") {return sf::Keyboard::K;}
    else if (s=="L") {return sf::Keyboard::L;}
    else if (s=="M") {return sf::Keyboard::M;}
    else if (s=="N") {return sf::Keyboard::N;}
    else if (s=="O") {return sf::Keyboard::O;}
    else if (s=="P") {return sf::Keyboard::P;}
    else if (s=="Q") {return sf::Keyboard::Q;}
    else if (s=="R") {return sf::Keyboard::R;}
    else if (s=="S") {return sf::Keyboard::S;}
    else if (s=="T") {return sf::Keyboard::T;}
    else if (s=="U") {return sf::Keyboard::U;}
    else if (s=="V") {return sf::Keyboard::V;}
    else if (s=="W") {return sf::Keyboard::W;}
    else if (s=="X") {return sf::Keyboard::X;}
    else if (s=="Y") {return sf::Keyboard::Y;}
    else if (s=="Z") {return sf::Keyboard::Z;}
    else if (s=="Num0") {return sf::Keyboard::Num0;}
    else if (s=="Num1") {return sf::Keyboard::Num1;}
    else if (s=="Num2") {return sf::Keyboard::Num2;}
    else if (s=="Num3") {return sf::Keyboard::Num3;}
    else if (s=="Num4") {return sf::Keyboard::Num4;}
    else if (s=="Num5") {return sf::Keyboard::Num5;}
    else if (s=="Num6") {return sf::Keyboard::Num6;}
    else if (s=="Num7") {return sf::Keyboard::Num7;}
    else if (s=="Num8") {return sf::Keyboard::Num8;}
    else if (s=="Num9") {return sf::Keyboard::Num9;}
    else if (s=="Escape") {return sf::Keyboard::Escape;}
    else if (s=="LControl") {return sf::Keyboard::LControl;}
    else if (s=="LShift") {return sf::Keyboard::LShift;}
    else if (s=="LAlt") {return sf::Keyboard::LAlt;}
    else if (s=="LSystem") {return sf::Keyboard::LSystem;}
    else if (s=="RControl") {return sf::Keyboard::RControl;}
    else if (s=="RShift") {return sf::Keyboard::RShift;}
    else if (s=="RAlt") {return sf::Keyboard::RAlt;}
    else if (s=="RSystem") {return sf::Keyboard::RSystem;}
    else if (s=="Menu") {return sf::Keyboard::Menu;}
    else if (s=="LBracket") {return sf::Keyboard::LBracket;}
    else if (s=="RBracket") {return sf::Keyboard::RBracket;}
    else if (s=="Semicolon") {return sf::Keyboard::Semicolon;}
    else if (s=="Comma") {return sf::Keyboard::Comma;}
    else if (s=="Period") {return sf::Keyboard::Period;}
    else if (s=="Quote") {return sf::Keyboard::Quote;}
    else if (s=="Slash") {return sf::Keyboard::Slash;}
    else if (s=="Backslash") {return sf::Keyboard::Backslash;}
    else if (s=="Tilde") {return sf::Keyboard::Tilde;}
    else if (s=="Equal") {return sf::Keyboard::Equal;}
    else if (s=="Hyphen") {return sf::Keyboard::Hyphen;}
    else if (s=="Space") {return sf::Keyboard::Space;}
    else if (s=="Enter") {return sf::Keyboard::Enter;}
    else if (s=="Backspace") {return sf::Keyboard::Backspace;}
    else if (s=="Tab") {return sf::Keyboard::Tab;}
    else if (s=="PageUp") {return sf::Keyboard::PageUp;}
    else if (s=="PageDown") {return sf::Keyboard::PageDown;}
    else if (s=="End") {return sf::Keyboard::End;}
    else if (s=="Home") {return sf::Keyboard::Home;}
    else if (s=="Insert") {return sf::Keyboard::Insert;}
    else if (s=="Delete") {return sf::Keyboard::Delete;}
    else if (s=="Add") {return sf::Keyboard::Add;}
    else if (s=="Subtract") {return sf::Keyboard::Subtract;}
    else if (s=="Multiply") {return sf::Keyboard::Multiply;}
    else if (s=="Divide") {return sf::Keyboard::Divide;}
    else if (s=="Left") {return sf::Keyboard::Left;}
    else if (s=="Right") {return sf::Keyboard::Right;}
    else if (s=="Up") {return sf::Keyboard::Up;}
    else if (s=="Down") {return sf::Keyboard::Down;}
    else if (s=="Numpad0") {return sf::Keyboard::Numpad0;}
    else if (s=="Numpad1") {return sf::Keyboard::Numpad1;}
    else if (s=="Numpad2") {return sf::Keyboard::Numpad2;}
    else if (s=="Numpad3") {return sf::Keyboard::Numpad3;}
    else if (s=="Numpad4") {return sf::Keyboard::Numpad4;}
    else if (s=="Numpad5") {return sf::Keyboard::Numpad5;}
    else if (s=="Numpad6") {return sf::Keyboard::Numpad6;}
    else if (s=="Numpad7") {return sf::Keyboard::Numpad7;}
    else if (s=="Numpad8") {return sf::Keyboard::Numpad8;}
    else if (s=="Numpad9") {return sf::Keyboard::Numpad9;}
    else if (s=="F1") {return sf::Keyboard::F1;}
    else if (s=="F2") {return sf::Keyboard::F2;}
    else if (s=="F3") {return sf::Keyboard::F3;}
    else if (s=="F4") {return sf::Keyboard::F4;}
    else if (s=="F5") {return sf::Keyboard::F5;}
    else if (s=="F6") {return sf::Keyboard::F6;}
    else if (s=="F7") {return sf::Keyboard::F7;}
    else if (s=="F8") {return sf::Keyboard::F8;}
    else if (s=="F9") {return sf::Keyboard::F9;}
    else if (s=="F10") {return sf::Keyboard::F10;}
    else if (s=="F11") {return sf::Keyboard::F11;}
    else if (s=="F12") {return sf::Keyboard::F12;}
    else if (s=="F13") {return sf::Keyboard::F13;}
    else if (s=="F14") {return sf::Keyboard::F14;}
    else if (s=="F15") {return sf::Keyboard::F15;}
    else if (s=="Pause") {return sf::Keyboard::Pause;}
    else if (s=="KeyCount") {return sf::Keyboard::KeyCount;}
    else if (s=="Unknown") {return sf::Keyboard::Unknown;}
}

#endif // USER-CONTROLS_H_INCLUDED
