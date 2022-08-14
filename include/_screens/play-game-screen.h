#ifndef PLAY-GAME-SCREEN_H_INCLUDED
#define PLAY-GAME-SCREEN_H_INCLUDED

class GameScreen : public GameState
{
    private:
        sf::Font _thinfont;
        sf::Font _boldfont;
        int gametype=0;
        //gametype 0 - multiplayer host.
        //gametype 1 - multiplayer but not host.
        //gametype 2 - singleplayer vs. AI.
        //gametype 3 - singleplayer lineup.
        Server server=Server();
        std::string targetip;
        unsigned short targetport;
        sf::TcpSocket socket;
        sf::Packet packet;
        sf::Uint16 packetId=0;
        sf::Uint32 resultsize=0;
        std::array<double,66> temp;

        sf::Font _font;
        sf::Font _scorefont;
        std::string p1name="PLAYER 1";
        std::string p2name="PLAYER 2";
        int p1score=0;
        int p2score=0;
        int p1frames=0;
        int p2frames=0;
        int framesbestof=35;
        bool placing_white=true;
        bool touching=false;
        bool isyourturn=true;
        bool isfoul=false;
        bool ismiss=false;
        bool isredon=true;
        int nom_colour_order=0; //1 yellow,2 green, etc.
        bool isfreeball=false;

        sf::Text ballontitle;

        bool change=true;
        bool done=true;
        int t=0;
        double deltat=0.;
        bool typing=false;

        double power=50.;
        double dist;
        double angle;

        bool ispaused=false;
        double pausea;
        bool ispausepressed=false;
        sf::RectangleShape pauserect;
        sf::Text pausetext;

        bool gameover=false;
        sf::RectangleShape gameoverrect;
        sf::Text gameovertext;

        //stats.
        sf::Text stats_title;
        std::vector<sf::Text> stats_text;
        std::vector<double> stats_text_y;
        int p1_highbreak=0;
        int p2_highbreak=0;
        int p1_centuries=0;
        int p2_centuries=0;

        sf::Vector2f spin_selector_pos;
        sf::Vector2f spin_dot_pos;

        sf::RectangleShape panel_background;
        sf::CircleShape spin_selector;
        sf::CircleShape spin_dot;
        sf::RectangleShape power_bar[100]={};
        sf::RectangleShape power_outline;
        sf::CircleShape elevation_border;
        sf::CircleShape elevation_ball;
        sf::RectangleShape elevation_pointer;
        sf::Text elevation_display;
        sf::Text textframesbestof;
        sf::Text textp1frames;
        sf::Text textp2frames;
        sf::Text textp1score;
        sf::Text textp2score;
        sf::Text textp1name;
        sf::Text textp2name;
        sf::CircleShape p1pointer;
        sf::CircleShape p2pointer;
        sf::RectangleShape framescoresrect;
        sf::RectangleShape p1scorerect;
        sf::RectangleShape p2scorerect;
        sf::RectangleShape p1namerect;
        sf::RectangleShape p2namerect;

        RedArrow redarrow;
        Cue cue;
        Ball balls[22];
        Cushion cushions[6];
        sf::CircleShape spots[6];
        sf::VertexArray baulkline;
        sf::VertexArray baulkcircle;

        Eigen::MatrixXd test=Eigen::MatrixXd::Constant(46,1,0.0);
        std::vector<std::array<double,66> > result;
        std::array<std::vector<double>,3> predict;

        sf::VertexArray cuetraj;
        sf::VertexArray cuetraj2;
        sf::VertexArray obtraj;
        sf::CircleShape ghostball;
        sf::CircleShape ghostball2;

        NominateBall nomballs[6];
        sf::RectangleShape nomback;
        sf::Text nomtext;

        static const int numlines=4;
        std::array<std::string,numlines> logstrings={};
        std::array<sf::Text,numlines> logtext;
        sf::RectangleShape logback;

        double ds=0.;
        sf::Vector2f pos;

    public:
        //set up the objects.
        GameScreen(double sfac=dfactor,int kind=0,std::string ip="", unsigned short port=50000,std::string name="PLAYER 1") : GameState(sfac)
        {
            gametype=kind;

            p1name=name;

            targetip=ip;
            targetport=port;

            if (gametype<2)
            {
                //multiplayer.
                if (socket.connect(sf::IpAddress(targetip),targetport)!=sf::Socket::Done)
                {
                    //error.
                    std::cout << "Could not connect to self server" << std::endl;
                }
                socket.setBlocking(false);
            }
            else
            {
                if (gametype==2)
                {
                    p2name="AI";
                }
                else
                {
                    p2name="N/A";
                    framesbestof=1;
                }
            }

            //set up some stuff depending on the gametype here.

            _background=sf::Color(baizecolour[0],baizecolour[1],baizecolour[2]);

            if (!_thinfont.loadFromFile(_thinFontFile))
            {
                std::cout << "Error loading font." << std::endl;
            }
            if (!_boldfont.loadFromFile(_boldFontFile))
            {
                std::cout << "Error loading font." << std::endl;
            }

            sf::FloatRect textrect;

            //set up pause rect.
            pauserect.setSize(sf::Vector2f(_sfac*raw_width,_sfac*raw_height));
            pauserect.setPosition(sf::Vector2f(0.,-_sfac*raw_height));
            pauserect.setFillColor(sf::Color(100,100,100,150));

            pausetext.setFont(_thinfont);
            pausetext.setCharacterSize(int(_sfac*raw_height*0.1));
            pausetext.setString("Paused");
            textrect=pausetext.getLocalBounds();
            pausetext.setOrigin(sf::Vector2f(textrect.left+0.5*textrect.width,textrect.top+0.5*textrect.height));
            pausetext.setPosition(sf::Vector2f(_sfac*raw_width*0.5,_sfac*raw_height*0.15-_sfac*raw_height));
            pausetext.setFillColor(sf::Color(255,255,255));

            //set up game over screen.
            gameoverrect.setSize(sf::Vector2f(_sfac*raw_width,_sfac*raw_height));
            gameoverrect.setPosition(sf::Vector2f(0.,-_sfac*raw_height));
            gameoverrect.setFillColor(sf::Color(100,100,100,150));

            gameovertext.setFont(_thinfont);
            gameovertext.setCharacterSize(int(_sfac*raw_height*0.1));
            gameovertext.setString("Game over");
            textrect=gameovertext.getLocalBounds();
            gameovertext.setOrigin(sf::Vector2f(textrect.left+0.5*textrect.width,textrect.top+0.5*textrect.height));
            gameovertext.setPosition(sf::Vector2f(_sfac*raw_width*0.5,_sfac*raw_height*0.15-_sfac*raw_height));
            gameovertext.setFillColor(sf::Color(255,255,255));

            int sh=int(_sfac*raw_height*0.04);

            stats_title.setFont(_thinfont);
            stats_title.setCharacterSize(int(sh*1.2));
            stats_title.setString("Overall Stats");
            textrect=stats_title.getLocalBounds();
            stats_title.setOrigin(sf::Vector2f(textrect.left+0.5*textrect.width,textrect.top+0.5*textrect.height));
            stats_title.setPosition(sf::Vector2f(_sfac*raw_width*0.5,_sfac*raw_height*0.25-_sfac*raw_height));
            stats_title.setFillColor(sf::Color(255,255,255));

            stats_text.push_back(sf::Text());
            stats_text.push_back(sf::Text());
            stats_text.push_back(sf::Text());
            stats_text.push_back(sf::Text());
            stats_text.push_back(sf::Text());
            stats_text.push_back(sf::Text());
            stats_text.push_back(sf::Text());
            stats_text.push_back(sf::Text());
            stats_text.push_back(sf::Text());
            stats_text.push_back(sf::Text());
            stats_text.push_back(sf::Text());
            stats_text.push_back(sf::Text());

            double starth=0.35*_sfac*raw_height;
            double soffset=0.15*_sfac*raw_width;

            for (int i=0;i<stats_text.size();i++)
            {
                stats_text[i].setFont(_thinfont);
                stats_text[i].setCharacterSize(sh);

                if (i==0) {stats_text[i].setString("versus");}
                else if (i==1) {stats_text[i].setString(p1name);}
                else if (i==2) {stats_text[i].setString(p2name);}
                else if (i==3) {stats_text[i].setString("Frames");}
                else if (i==4) {stats_text[i].setString(std::to_string(p1frames));}
                else if (i==5) {stats_text[i].setString(std::to_string(p2frames));}
                else if (i==6) {stats_text[i].setString("High break");}
                else if (i==7) {stats_text[i].setString(std::to_string(p1_highbreak));}
                else if (i==8) {stats_text[i].setString(std::to_string(p2_highbreak));}
                else if (i==9) {stats_text[i].setString("Centuries");}
                else if (i==10) {stats_text[i].setString(std::to_string(p1_centuries));}
                else if (i==11) {stats_text[i].setString(std::to_string(p2_centuries));}

                textrect=stats_text[i].getLocalBounds();
                stats_text[i].setOrigin(sf::Vector2f(textrect.left+0.5*textrect.width,textrect.top+0.5*textrect.height));
                if (i%3==0)
                {
                    stats_text[i].setPosition(sf::Vector2f(0.5*_sfac*raw_width,-_sfac*raw_height+starth+((i/3)*1.4)*sh));
                }
                else if (i%3==1)
                {
                    stats_text[i].setPosition(sf::Vector2f(0.5*_sfac*raw_width-soffset,-_sfac*raw_height+starth+((i/3)*1.4)*sh));
                }
                else if (i%3==2)
                {
                    stats_text[i].setPosition(sf::Vector2f(0.5*_sfac*raw_width+soffset,-_sfac*raw_height+starth+((i/3)*1.4)*sh));
                }
                stats_text[i].setFillColor(sf::Color(255,255,255));

                stats_text_y.push_back(stats_text[i].getPosition().y);
            }

            pausea=2.*_sfac*raw_height;

            _buttons.push_back(RectButton());
            _buttons.push_back(RectButton());
            _buttons.push_back(RectButton());
            _buttons.push_back(RectButton());

            double buttonwidth=_sfac*raw_width*0.15;
            int textcolour1[4]={255,255,255,255};
            int textcolour2[4]={255,255,255,255};
            int textcolour3[4]={255,255,255,255};
            int outlinecolour1[4]={200,200,200,255};
            int outlinecolour2[4]={255,255,255,255};
            int outlinecolour3[4]={255,255,255,255};
            int colour1[4]={100,100,100,150};
            int colour2[4]={100,100,100,150};
            int colour3[4]={169,169,169,200};

            for (int i=0;i<2;i++)
            {
                _buttons[i]._shape.setSize(sf::Vector2f(buttonwidth,buttonwidth/_buttons[i]._ratio));
                _buttons[i]._shape.setOrigin(0.5*buttonwidth,0.5*buttonwidth/_buttons[i]._ratio);
                _buttons[i]._shape.setOutlineThickness(_buttons[i]._absoutlinethickness);
                if (i==0)
                {
                    _buttons[i]._shape.setPosition(sf::Vector2f(0.5*_sfac*raw_width,-_sfac*raw_height+0.50*_sfac*raw_height+(i*1.2)*buttonwidth/_buttons[i]._ratio));
                }
                else if (i==1)
                {
                    _buttons[i]._shape.setPosition(sf::Vector2f(0.5*_sfac*raw_width,-_sfac*raw_height+0.75*_sfac*raw_height));
                }

                if (!_buttons[i]._font.loadFromFile(_thinFontFile)) {std::cout << "Error loading font." << std::endl;}
                _buttons[i]._text.setFont(_buttons[i]._font);
                _buttons[i]._text.setCharacterSize(int(buttonwidth*_buttons[i]._textfactor/_buttons[i]._ratio));
                if (i==0)
                {
                    _buttons[i]._text.setString("Exit");
                    _buttons[i]._target="Quit";
                }
                if (i==1)
                {
                    _buttons[i]._text.setString("Exit");
                    _buttons[i]._target="Quit";
                }
                textrect=_buttons[i]._text.getLocalBounds();
                _buttons[i]._text.setOrigin(sf::Vector2f(textrect.left+0.5*textrect.width,textrect.top+0.5*textrect.height));

                if (i==0)
                {
                    _buttons[i]._text.setPosition(sf::Vector2f(0.5*_sfac*raw_width,-_sfac*raw_height+0.50*_sfac*raw_height+(i*1.2)*buttonwidth/_buttons[i]._ratio));
                }
                else if (i==1)
                {
                    _buttons[i]._text.setPosition(sf::Vector2f(0.5*_sfac*raw_width,-_sfac*raw_height+0.75*_sfac*raw_height));
                }

                _buttons[i]._colour1=sf::Color(colour1[0],colour1[1],colour1[2],colour1[3]);
                _buttons[i]._colour2=sf::Color(colour2[0],colour2[1],colour2[2],colour2[3]);
                _buttons[i]._colour3=sf::Color(colour3[0],colour3[1],colour3[2],colour3[3]);

                _buttons[i]._textcolour1=sf::Color(textcolour1[0],textcolour1[1],textcolour1[2],textcolour1[3]);
                _buttons[i]._textcolour2=sf::Color(textcolour2[0],textcolour2[1],textcolour2[2],textcolour2[3]);
                _buttons[i]._textcolour3=sf::Color(textcolour3[0],textcolour3[1],textcolour3[2],textcolour3[3]);

                _buttons[i]._outlinecolour1=sf::Color(outlinecolour1[0],outlinecolour1[1],outlinecolour1[2],outlinecolour1[3]);
                _buttons[i]._outlinecolour2=sf::Color(outlinecolour2[0],outlinecolour2[1],outlinecolour2[2],outlinecolour2[3]);
                _buttons[i]._outlinecolour3=sf::Color(outlinecolour3[0],outlinecolour3[1],outlinecolour3[2],outlinecolour3[3]);

                _buttons[i]._shape.setFillColor(_buttons[i]._colour1);
                _buttons[i]._shape.setOutlineColor(_buttons[i]._outlinecolour1);
                _buttons[i]._text.setFillColor(_buttons[i]._textcolour1);
            }

            double sbwidth=(_sfac*raw_height*panel_ratio/(1.+panel_ratio))*6.;
            double sbheight=(_sfac*raw_height*panel_ratio/(1.+panel_ratio))*0.2;

            buttonwidth=0.12*_sfac*raw_width;
            for (int i=2;i<4;i++)
            {
                _buttons[i]._shape.setSize(sf::Vector2f(buttonwidth,sbheight));
                _buttons[i]._shape.setOrigin(sf::Vector2f(0.5*buttonwidth,0.5*sbheight));
                _buttons[i]._shape.setOutlineThickness(_buttons[i]._absoutlinethickness);
                _buttons[i]._shape.setPosition(sf::Vector2f(0.5*_sfac*raw_width+0.5*sbwidth-0.5*buttonwidth,0.92*_sfac*raw_height+1.2*(i-2)*sbheight));

                if (!_buttons[i]._font.loadFromFile(_thinFontFile)) {std::cout << "Error loading font." << std::endl;}
                _buttons[i]._text.setFont(_buttons[i]._font);
                _buttons[i]._text.setCharacterSize(int(0.7*sbheight));
                if (i==2)
                {
                    _buttons[i]._text.setString("Concede frame");
                    _buttons[i]._target="Concedeframe";
                }
                else if (i==3)
                {
                    _buttons[i]._text.setString("Concede match");
                    _buttons[i]._target="Concedematch";
                }
                textrect=_buttons[i]._text.getLocalBounds();
                _buttons[i]._text.setOrigin(sf::Vector2f(textrect.left+0.5*textrect.width,textrect.top+0.5*textrect.height));
                _buttons[i]._text.setPosition(sf::Vector2f(0.5*_sfac*raw_width+0.5*sbwidth-0.5*buttonwidth,0.92*_sfac*raw_height+1.2*(i-2)*sbheight));

                _buttons[i]._colour1=sf::Color(colour1[0],colour1[1],colour1[2],colour1[3]);
                _buttons[i]._colour2=sf::Color(colour2[0],colour2[1],colour2[2],colour2[3]);
                _buttons[i]._colour3=sf::Color(colour3[0],colour3[1],colour3[2],colour3[3]);

                _buttons[i]._textcolour1=sf::Color(textcolour1[0],textcolour1[1],textcolour1[2],textcolour1[3]);
                _buttons[i]._textcolour2=sf::Color(textcolour2[0],textcolour2[1],textcolour2[2],textcolour2[3]);
                _buttons[i]._textcolour3=sf::Color(textcolour3[0],textcolour3[1],textcolour3[2],textcolour3[3]);

                _buttons[i]._outlinecolour1=sf::Color(outlinecolour1[0],outlinecolour1[1],outlinecolour1[2],outlinecolour1[3]);
                _buttons[i]._outlinecolour2=sf::Color(outlinecolour2[0],outlinecolour2[1],outlinecolour2[2],outlinecolour2[3]);
                _buttons[i]._outlinecolour3=sf::Color(outlinecolour3[0],outlinecolour3[1],outlinecolour3[2],outlinecolour3[3]);

                _buttons[i]._shape.setFillColor(_buttons[i]._colour1);
                _buttons[i]._shape.setOutlineColor(_buttons[i]._outlinecolour1);
                _buttons[i]._text.setFillColor(_buttons[i]._textcolour1);
            }

            //set up panel.
            panel_background=sf::RectangleShape(sf::Vector2f((_sfac*raw_width),(_sfac*raw_height*panel_ratio/(1.+panel_ratio))));
            panel_background.setPosition(sf::Vector2f(0.,(_sfac*raw_height/(1.+panel_ratio))));
            panel_background.setFillColor(sf::Color(0,0,0));

            spin_selector.setRadius((_sfac*raw_height*panel_ratio/(1.+panel_ratio))*0.4);
            spin_selector.setOrigin((_sfac*raw_height*panel_ratio/(1.+panel_ratio))*0.4,(_sfac*raw_height*panel_ratio/(1.+panel_ratio))*0.4);
            spin_selector.setPosition(sf::Vector2f(0.5*(_sfac*raw_height*panel_ratio/(1.+panel_ratio)),(_sfac*raw_height/(1.+panel_ratio))+0.5*(_sfac*raw_height*panel_ratio/(1.+panel_ratio))));
            spin_selector.setFillColor(sf::Color(255,255,255));
            spin_selector.setPointCount(200);

            spin_dot.setRadius((_sfac*raw_height*panel_ratio/(1.+panel_ratio))*0.4*0.005/0.02625);
            spin_dot.setOrigin(sf::Vector2f((_sfac*raw_height*panel_ratio/(1.+panel_ratio))*0.4*0.005/0.02625,(_sfac*raw_height*panel_ratio/(1.+panel_ratio))*0.4*0.005/0.02625));
            spin_dot.setPosition(sf::Vector2f(0.5*(_sfac*raw_height*panel_ratio/(1.+panel_ratio)),(_sfac*raw_height/(1.+panel_ratio))+0.5*(_sfac*raw_height*panel_ratio/(1.+panel_ratio))));
            spin_dot.setFillColor(sf::Color(255,0,0));
            spin_dot.setPointCount(100);

            for (int i=0;i<100;i++)
            {
                power_bar[i]=sf::RectangleShape(sf::Vector2f(0.1*(_sfac*raw_height*panel_ratio/(1.+panel_ratio)),0.008*(_sfac*raw_height*panel_ratio/(1.+panel_ratio))));
                power_bar[i].setPosition(sf::Vector2f(1.1*(_sfac*raw_height*panel_ratio/(1.+panel_ratio)),(_sfac*raw_height/(1.+panel_ratio))+0.1*(_sfac*raw_height*panel_ratio/(1.+panel_ratio))+i*0.008*(_sfac*raw_height*panel_ratio/(1.+panel_ratio))));
                power_bar[i].setFillColor(sf::Color(255,int(floor(2.55*i)),0));
            }
            sf::Color c;
            for (int i=0;i<100;i++)
            {
                c=power_bar[99-i].getFillColor();
                if (i<int(floor(power))) {power_bar[99-i].setFillColor(sf::Color(c.r,c.g,c.b,255));}
                else {power_bar[99-i].setFillColor(sf::Color(c.r,c.g,c.b,0));}
            }
            double outline=0.01*(_sfac*raw_height*panel_ratio/(1.+panel_ratio));
            power_outline=sf::RectangleShape(sf::Vector2f(0.1*(_sfac*raw_height*panel_ratio/(1.+panel_ratio))+2*outline,0.8*(_sfac*raw_height*panel_ratio/(1.+panel_ratio))+2*outline));
            power_outline.setPosition(sf::Vector2f(1.1*(_sfac*raw_height*panel_ratio/(1.+panel_ratio))-outline,(_sfac*raw_height/(1.+panel_ratio))+0.1*(_sfac*raw_height*panel_ratio/(1.+panel_ratio))-outline));
            power_outline.setFillColor(sf::Color(0,0,0));
            power_outline.setOutlineThickness(outline);
            power_outline.setOutlineColor(sf::Color(255,255,255));

            elevation_border.setRadius((_sfac*raw_height*panel_ratio/(1.+panel_ratio))*0.4);
            elevation_border.setOrigin((_sfac*raw_height*panel_ratio/(1.+panel_ratio))*0.4,(_sfac*raw_height*panel_ratio/(1.+panel_ratio))*0.4);
            elevation_border.setPosition(sf::Vector2f((_sfac*raw_width)-0.5*(_sfac*raw_height*panel_ratio/(1.+panel_ratio)),(_sfac*raw_height/(1.+panel_ratio))+0.5*(_sfac*raw_height*panel_ratio/(1.+panel_ratio))));
            elevation_border.setFillColor(sf::Color(0,0,0));
            elevation_border.setOutlineThickness(outline);
            elevation_border.setOutlineColor(sf::Color(255,255,255));
            elevation_border.setPointCount(200);

            elevation_ball.setRadius(0.2*(_sfac*raw_height*panel_ratio/(1.+panel_ratio)));
            elevation_ball.setOrigin(0.2*(_sfac*raw_height*panel_ratio/(1.+panel_ratio)),0.2*(_sfac*raw_height*panel_ratio/(1.+panel_ratio)));
            elevation_ball.setPosition(sf::Vector2f((_sfac*raw_width)-0.5*(_sfac*raw_height*panel_ratio/(1.+panel_ratio)),(_sfac*raw_height/(1.+panel_ratio))+0.5*(_sfac*raw_height*panel_ratio/(1.+panel_ratio))));
            elevation_ball.setFillColor(sf::Color(255,255,255));
            elevation_ball.setPointCount(150);

            elevation_pointer=sf::RectangleShape(sf::Vector2f(0.2*(_sfac*raw_height*panel_ratio/(1.+panel_ratio)),0.2*(_sfac*raw_height*panel_ratio/(1.+panel_ratio))*0.005/0.02625));
            elevation_pointer.setOrigin(0,0.1*(_sfac*raw_height*panel_ratio/(1.+panel_ratio))*0.005/0.02625);
            elevation_pointer.setPosition(sf::Vector2f((_sfac*raw_width)-0.3*(_sfac*raw_height*panel_ratio/(1.+panel_ratio)),(_sfac*raw_height/(1.+panel_ratio))+0.5*(_sfac*raw_height*panel_ratio/(1.+panel_ratio))));
            elevation_pointer.setFillColor(sf::Color(255,0,0));

            if (!_font.loadFromFile(_thinFontFile))
            {
                std::cout << "ERROR. Couldn't load font!" << std::endl;
            }
            elevation_display.setFont(_font);
            elevation_display.setString("00°");
            elevation_display.setCharacterSize(int(0.15*(_sfac*raw_height*panel_ratio/(1.+panel_ratio))));
            elevation_display.setFillColor(sf::Color(0,0,0));
            sf::FloatRect textRect=elevation_display.getLocalBounds();
            elevation_display.setOrigin(sf::Vector2f(textRect.left+textRect.width/2.,textRect.top+textRect.height/2.));
            elevation_display.setPosition(sf::Vector2f((_sfac*raw_width)-0.5*(_sfac*raw_height*panel_ratio/(1.+panel_ratio)),(_sfac*raw_height/(1.+panel_ratio))+0.5*(_sfac*raw_height*panel_ratio/(1.+panel_ratio))));

            if (!_scorefont.loadFromFile(_boldFontFile))
            {
                std::cout << "ERROR. Couldn't load font!" << std::endl;
            }
            textframesbestof.setFont(_scorefont);
            textframesbestof.setString("("+std::to_string(framesbestof)+")");
            textframesbestof.setCharacterSize(int(sbheight*0.7));
            textframesbestof.setFillColor(sf::Color(255,255,255));
            textRect=textframesbestof.getLocalBounds();
            textframesbestof.setOrigin(sf::Vector2f(textRect.left+textRect.width/2.,textRect.top+textRect.height/2.));
            textframesbestof.setPosition(sf::Vector2f((_sfac*raw_width)*0.5,(_sfac*raw_height/(1.+panel_ratio))+sbheight));

            textp1frames.setFont(_scorefont);
            textp1frames.setString(std::to_string(p1frames));
            textp1frames.setCharacterSize(int(sbheight*0.7));
            textp1frames.setFillColor(sf::Color(255,255,255));
            textRect=textp1frames.getLocalBounds();
            textp1frames.setOrigin(sf::Vector2f(textRect.left+textRect.width/2.,textRect.top+textRect.height/2.));
            textp1frames.setPosition(sf::Vector2f((_sfac*raw_width)*0.5-0.06*sbwidth,(_sfac*raw_height/(1.+panel_ratio))+sbheight));

            textp2frames.setFont(_scorefont);
            textp2frames.setString(std::to_string(p2frames));
            textp2frames.setCharacterSize(int(sbheight*0.7));
            textp2frames.setFillColor(sf::Color(255,255,255));
            textRect=textp2frames.getLocalBounds();
            textp2frames.setOrigin(sf::Vector2f(textRect.left+textRect.width/2.,textRect.top+textRect.height/2.));
            textp2frames.setPosition(sf::Vector2f((_sfac*raw_width)*0.5+0.06*sbwidth,(_sfac*raw_height/(1.+panel_ratio))+sbheight));

            textp1score.setFont(_scorefont);
            textp1score.setString(std::to_string(p1score));
            textp1score.setCharacterSize(int(sbheight*0.7));
            textp1score.setFillColor(sf::Color(0,0,0));
            textRect=textp1score.getLocalBounds();
            textp1score.setOrigin(sf::Vector2f(textRect.left+textRect.width/2.,textRect.top+textRect.height/2.));
            textp1score.setPosition(sf::Vector2f((_sfac*raw_width)*0.5-0.1*sbwidth-0.5*0.15*sbwidth,(_sfac*raw_height/(1.+panel_ratio))+sbheight));

            textp2score.setFont(_scorefont);
            textp2score.setString(std::to_string(p2score));
            textp2score.setCharacterSize(int(sbheight*0.7));
            textp2score.setFillColor(sf::Color(0,0,0));
            textRect=textp2score.getLocalBounds();
            textp2score.setOrigin(sf::Vector2f(textRect.left+textRect.width/2.,textRect.top+textRect.height/2.));
            textp2score.setPosition(sf::Vector2f((_sfac*raw_width)*0.5+0.1*sbwidth+0.5*0.15*sbwidth,(_sfac*raw_height/(1.+panel_ratio))+sbheight));

            textp1name.setFont(_scorefont);
            textp1name.setString(p1name);
            textp1name.setCharacterSize(int(sbheight*0.7));
            textp1name.setFillColor(sf::Color(0,0,0));
            textRect=textp1name.getLocalBounds();
            textp1name.setOrigin(sf::Vector2f(0.,textRect.top+textRect.height/2.));
            textp1name.setPosition(sf::Vector2f((_sfac*raw_width)*0.5-0.47*sbwidth,(_sfac*raw_height/(1.+panel_ratio))+sbheight));

            textp2name.setFont(_scorefont);
            textp2name.setString(p2name);
            textp2name.setCharacterSize(int(sbheight*0.7));
            textp2name.setFillColor(sf::Color(0,0,0));
            textRect=textp2name.getLocalBounds();
            textp2name.setOrigin(sf::Vector2f(textRect.left+textRect.width,textRect.top+textRect.height/2.));
            textp2name.setPosition(sf::Vector2f((_sfac*raw_width)*0.5+0.47*sbwidth,(_sfac*raw_height/(1.+panel_ratio))+sbheight));

            p1pointer.setRadius(sbheight*0.2);
            p1pointer.setPointCount(3);
            p1pointer.setFillColor(sf::Color(0,0,0));
            p1pointer.setOrigin(sbheight*0.2,sbheight*0.2);
            p1pointer.setRotation(-90.);
            p1pointer.setPosition(sf::Vector2f((_sfac*raw_width)*0.5-0.105*sbwidth,(_sfac*raw_height/(1.+panel_ratio))+sbheight));

            p2pointer.setRadius(sbheight*0.2);
            p2pointer.setPointCount(3);
            p2pointer.setFillColor(sf::Color(0,0,0,0));
            p2pointer.setOrigin(sbheight*0.2,sbheight*0.2);
            p2pointer.setRotation(90.);
            p2pointer.setPosition(sf::Vector2f((_sfac*raw_width)*0.5+0.105*sbwidth,(_sfac*raw_height/(1.+panel_ratio))+sbheight));

            framescoresrect=sf::RectangleShape(sf::Vector2f(sbwidth*0.2,sbheight));
            framescoresrect.setOrigin(sbwidth*0.1,0.5*sbheight);
            framescoresrect.setPosition(sf::Vector2f((_sfac*raw_width)*0.5,(_sfac*raw_height/(1.+panel_ratio))+sbheight));
            framescoresrect.setFillColor(sf::Color(51,153,255));

            p1scorerect=sf::RectangleShape(sf::Vector2f(sbwidth*0.15,sbheight));
            p1scorerect.setOrigin(sbwidth*0.15*0.5,0.5*sbheight);
            p1scorerect.setPosition(sf::Vector2f((_sfac*raw_width)*0.5-0.1*sbwidth-0.15*0.5*sbwidth,(_sfac*raw_height/(1.+panel_ratio))+sbheight));
            p1scorerect.setFillColor(sf::Color(255,255,255));

            p2scorerect=sf::RectangleShape(sf::Vector2f(sbwidth*0.15,sbheight));
            p2scorerect.setOrigin(sbwidth*0.15*0.5,0.5*sbheight);
            p2scorerect.setPosition(sf::Vector2f((_sfac*raw_width)*0.5+0.1*sbwidth+0.15*0.5*sbwidth,(_sfac*raw_height/(1.+panel_ratio))+sbheight));
            p2scorerect.setFillColor(sf::Color(255,255,255));

            p1namerect=sf::RectangleShape(sf::Vector2f(sbwidth*0.25,sbheight));
            p1namerect.setOrigin(sbwidth*0.25*0.5,0.5*sbheight);
            p1namerect.setPosition(sf::Vector2f((_sfac*raw_width)*0.5-0.25*sbwidth-0.25*0.5*sbwidth,(_sfac*raw_height/(1.+panel_ratio))+sbheight));
            p1namerect.setFillColor(sf::Color(236,228,0));

            p2namerect=sf::RectangleShape(sf::Vector2f(sbwidth*0.25,sbheight));
            p2namerect.setOrigin(sbwidth*0.25*0.5,0.5*sbheight);
            p2namerect.setPosition(sf::Vector2f((_sfac*raw_width)*0.5+0.25*sbwidth+0.25*0.5*sbwidth,(_sfac*raw_height/(1.+panel_ratio))+sbheight));
            p2namerect.setFillColor(sf::Color(236,228,0));

            _shapes.push_back(&panel_background);
            _shapes.push_back(&spin_selector);
            _shapes.push_back(&spin_dot);
            _shapes.push_back(&power_outline);
            for (int i=0;i<100;i++)
            {
                _shapes.push_back(&power_bar[i]);
            }
            _shapes.push_back(&elevation_border);
            _shapes.push_back(&elevation_ball);
            _shapes.push_back(&elevation_pointer);
            _shapes.push_back(&framescoresrect);
            _shapes.push_back(&p1scorerect);
            _shapes.push_back(&p2scorerect);
            _shapes.push_back(&p1namerect);
            _shapes.push_back(&p2namerect);
            _shapes.push_back(&p1pointer);
            _shapes.push_back(&p2pointer);
            _shapes.push_back(&elevation_display);
            _shapes.push_back(&textframesbestof);
            _shapes.push_back(&textp1frames);
            _shapes.push_back(&textp2frames);
            _shapes.push_back(&textp1score);
            _shapes.push_back(&textp2score);
            _shapes.push_back(&textp1name);
            _shapes.push_back(&textp2name);

            ballontitle.setFont(_thinfont);
            ballontitle.setString("Ball on - Red");
            ballontitle.setCharacterSize(int(sbheight*0.7));
            textrect=ballontitle.getLocalBounds();
            ballontitle.setOrigin(sf::Vector2f(textrect.left,textrect.top));
            ballontitle.setPosition(sf::Vector2f((_sfac*raw_width)*0.5-0.5*sbwidth,_sfac*raw_height*0.89));
            ballontitle.setFillColor(sf::Color(255,255,255));
            _shapes.push_back(&ballontitle);

            nomtext.setFont(_thinfont);
            nomtext.setString("Nominate colour");
            nomtext.setCharacterSize(int(sbheight*0.7));
            textrect=nomtext.getLocalBounds();
            nomtext.setOrigin(sf::Vector2f(textrect.left,textrect.top));
            nomtext.setPosition(sf::Vector2f((_sfac*raw_width)*0.5-0.5*sbwidth,_sfac*raw_height*0.89+1.1*int(sbheight*0.7)));
            nomtext.setFillColor(sf::Color(255,255,255,100));
            _shapes.push_back(&nomtext);

            nomback.setSize(sf::Vector2f(((6.)*1.2*2.)*_sfac*raw_height*0.015,1.3*_sfac*raw_height*2.*0.015));
            nomback.setPosition(sf::Vector2f((_sfac*raw_width)*0.5-0.5*sbwidth,_sfac*raw_height*0.96-1.3*_sfac*raw_height*0.015));
            nomback.setFillColor(sf::Color(colour1[0],colour1[1],colour1[2],colour1[3]));
            nomback.setOutlineThickness(_buttons[0]._absoutlinethickness);
            nomback.setOutlineColor(sf::Color(outlinecolour1[0],outlinecolour1[1],outlinecolour1[2],outlinecolour1[3]));
            _shapes.push_back(&nomback);
            for (int i=0;i<6;i++)
            {
                nomballs[i]._shape.setRadius(_sfac*raw_height*0.015);
                nomballs[i]._shape.setOrigin(sf::Vector2f(_sfac*raw_height*0.015,_sfac*raw_height*0.015));
                nomballs[i]._shape.setPosition(sf::Vector2f((_sfac*raw_width)*0.5-0.5*sbwidth+((i+0.5)*1.2*2.)*_sfac*raw_height*0.015,_sfac*raw_height*0.96));
                if (i==0)
                {
                    nomballs[i]._colour=yellow_col;
                }
                else if (i==1)
                {
                    nomballs[i]._colour=green_col;
                }
                else if (i==2)
                {
                    nomballs[i]._colour=brown_col;
                }
                else if (i==3)
                {
                    nomballs[i]._colour=blue_col;
                }
                else if (i==4)
                {
                    nomballs[i]._colour=pink_col;
                }
                else if (i==5)
                {
                    nomballs[i]._colour=black_col;
                }
                nomballs[i]._shape.setFillColor(nomballs[i]._colour);
                nomballs[i]._shape.setOutlineThickness(0.);
                nomballs[i]._shape.setOutlineColor(sf::Color(255,255,255,100));

                _shapes.push_back(&nomballs[i]._shape);
            }

            double logheight=(_sfac*raw_height*panel_ratio/(1.+panel_ratio))*0.6*0.75;
            logback.setSize(sf::Vector2f(0.5*sbwidth,logheight));
            logback.setPosition(sf::Vector2f(0.5*_sfac*raw_width-0.25*sbwidth,0.89*_sfac*raw_height));
            logback.setOutlineThickness(_buttons[0]._absoutlinethickness);
            logback.setFillColor(_buttons[0]._colour1);
            logback.setOutlineColor(_buttons[0]._outlinecolour1);
            _shapes.push_back(&logback);

            double cheight=logheight/(numlines*1.2);
            for (int i=0;i<numlines;i++)
            {
                logstrings[i]="";
                logtext[i].setFont(_thinfont);
                logtext[i].setCharacterSize(int(cheight));
                logtext[i].setString(logstrings[i]);
                textrect=logtext[i].getLocalBounds();
                logtext[i].setOrigin(sf::Vector2f(textrect.left,textrect.top));
                logtext[i].setPosition(sf::Vector2f(0.5*_sfac*raw_width-0.25*sbwidth+0.1*cheight,0.89*_sfac*raw_height+0.1*cheight+1.2*i*cheight));
                logtext[i].setFillColor(sf::Color(255,255,255));
                _shapes.push_back(&logtext[i]);
            }

            _inputboxes.push_back(InputBox());

            if (gametype<2)
            {
                for (int i=0;i<1;i++)
                {
                    _inputboxes[i]._shape.setSize(sf::Vector2f(0.5*sbwidth,cheight/_inputboxes[i]._textfactor));
                    _inputboxes[i]._shape.setOrigin(sf::Vector2f(0.25*sbwidth,0.5*cheight/_inputboxes[i]._textfactor));
                    _inputboxes[i]._shape.setOutlineThickness(_inputboxes[i]._absoutlinethickness);
                    _inputboxes[i]._shape.setPosition(sf::Vector2f(0.5*_sfac*raw_width,0.89*_sfac*raw_height+1.05*logheight+0.5*cheight/_inputboxes[i]._textfactor));

                    if (!_inputboxes[i]._font.loadFromFile(_thinFontFile)) {std::cout << "Error loading font." << std::endl;}
                    _inputboxes[i]._text.setFont(_inputboxes[i]._font);
                    _inputboxes[i]._text.setCharacterSize(int(cheight));
                    _inputboxes[i]._text.setFillColor(sf::Color(255,255,255));
                    _inputboxes[i]._text.setString("");
                    textrect=_inputboxes[i]._text.getLocalBounds();
                    _inputboxes[i]._text.setOrigin(sf::Vector2f(textrect.left,textrect.top));
                    _inputboxes[i]._text.setPosition(sf::Vector2f(0.5*_sfac*raw_width-0.25*sbwidth+2.*_inputboxes[i]._absoutlinethickness,0.89*_sfac*raw_height+1.05*logheight+(cheight/_inputboxes[i]._textfactor)*(0.5*(1.-_inputboxes[i]._textfactor))));

                    if (i==0)
                    {
                        _inputboxes[i]._backtext.setString("Enter message...");
                    }
                    _inputboxes[i]._backtext.setFont(_inputboxes[i]._font);
                    _inputboxes[i]._backtext.setCharacterSize(int(cheight));
                    textrect=_inputboxes[i]._backtext.getLocalBounds();
                    _inputboxes[i]._backtext.setOrigin(sf::Vector2f(textrect.left,textrect.top+textrect.height*0.5));
                    _inputboxes[i]._backtext.setFillColor(sf::Color(255,255,255,150));
                    _inputboxes[i]._backtext.setPosition(sf::Vector2f(0.5*_sfac*raw_width-0.25*sbwidth+2.*_inputboxes[i]._absoutlinethickness,0.89*_sfac*raw_height+1.05*logheight+0.5*(cheight/_inputboxes[i]._textfactor)));

                    _inputboxes[i]._abscursorthickness=0.5;
                    _inputboxes[i]._cursor.setSize(sf::Vector2f(_inputboxes[i]._abscursorthickness,cheight));
                    _inputboxes[i]._cursor.setOrigin(sf::Vector2f(0.5*_inputboxes[i]._abscursorthickness,0.5*cheight));
                    _inputboxes[i]._cursor.setPosition(sf::Vector2f(0.5*_sfac*raw_width-0.25*sbwidth+2.*_inputboxes[i]._absoutlinethickness,0.89*_sfac*raw_height+1.05*logheight+0.5*(cheight/_inputboxes[i]._textfactor)));

                    _inputboxes[i]._outlinecolour1=sf::Color(outlinecolour1[0],outlinecolour1[1],outlinecolour1[2],outlinecolour1[3]);
                    _inputboxes[i]._outlinecolour2=sf::Color(outlinecolour2[0],outlinecolour2[1],outlinecolour2[2],outlinecolour2[3]);

                    _inputboxes[i]._shape.setFillColor(sf::Color(colour1[0],colour1[1],colour1[2],colour1[3]));
                    _inputboxes[i]._shape.setOutlineColor(_inputboxes[i]._outlinecolour1);
                    _inputboxes[i]._cursor.setFillColor(sf::Color(255,255,255,0));
                }
            }

            cushions[0]=Cushion(mpockets[0][0],mpockets[0][1]-0.156,pi,1,0);
            cushions[1]=Cushion(mpockets[1][0],mpockets[1][1]+0.156,0.0,1,0);
            cushions[2]=Cushion(cpockets[0][0],cpockets[0][1],0.0,0,1);
            cushions[3]=Cushion(cpockets[1][0],cpockets[1][1],pi/2,0,0);
            cushions[4]=Cushion(cpockets[2][0],cpockets[2][1],pi*3/2,0,0);
            cushions[5]=Cushion(cpockets[3][0],cpockets[3][1],pi,0,1);

            for (int i=0;i<6;i++)
            {
                _shapes.push_back(&cushions[i]._pocketshape);
            }
            for (int i=0;i<6;i++)
            {
                _shapes.push_back(&cushions[i]._shape);
                _shapes.push_back(&cushions[i]._railshape);
                _shapes.push_back(&cushions[i]._p1shape);
                _shapes.push_back(&cushions[i]._p2shape);
            }

            for (int i=0;i<6;i++)
            {
                spots[i].setRadius(spot_r*_sfac);
                spots[i].setOrigin(spot_r*_sfac,spot_r*_sfac);
                spots[i].setFillColor(sf::Color(255,255,255));
            }

            for (int i=0;i<22;i++)
            {
                balls[i]._shape.setOrigin(ball_radius*_sfac,ball_radius*_sfac);
            }

            //cueball.
            balls[0]._shape.setFillColor(sf::Color(255,255,255));
            balls[0]._x=cueball_break_x;
            balls[0]._y=cueball_break_y;
            balls[0]._order=1;
            //yellow.
            balls[1]._shape.setFillColor(yellow_col);
            balls[1]._x=yellow_x;
            balls[1]._y=yellow_y;
            spots[0].setPosition(sf::Vector2f(yellow_x*_sfac,(_sfac*raw_height/(1.+panel_ratio))-yellow_y*_sfac));
            balls[1]._order=2;
            //green.
            balls[2]._shape.setFillColor(green_col);
            balls[2]._x=green_x;
            balls[2]._y=green_y;
            spots[1].setPosition(sf::Vector2f(green_x*_sfac,(_sfac*raw_height/(1.+panel_ratio))-green_y*_sfac));
            balls[2]._order=3;
            //brown.
            balls[3]._shape.setFillColor(brown_col);
            balls[3]._x=brown_x;
            balls[3]._y=brown_y;
            spots[2].setPosition(sf::Vector2f(brown_x*_sfac,(_sfac*raw_height/(1.+panel_ratio))-brown_y*_sfac));
            balls[3]._order=4;
            //blue.
            balls[4]._shape.setFillColor(blue_col);
            balls[4]._x=blue_x;
            balls[4]._y=blue_y;
            spots[3].setPosition(sf::Vector2f(blue_x*_sfac,(_sfac*raw_height/(1.+panel_ratio))-blue_y*_sfac));
            balls[4]._order=5;
            //pink.
            balls[5]._shape.setFillColor(pink_col);
            balls[5]._x=pink_x;
            balls[5]._y=pink_y;
            spots[4].setPosition(sf::Vector2f(pink_x*_sfac,(_sfac*raw_height/(1.+panel_ratio))-pink_y*_sfac));
            balls[5]._order=6;
            //black.
            balls[6]._shape.setFillColor(black_col);
            balls[6]._x=black_x;
            balls[6]._y=black_y;
            spots[5].setPosition(sf::Vector2f(black_x*_sfac,(_sfac*raw_height/(1.+panel_ratio))-black_y*_sfac));
            balls[6]._order=7;

            double x;
            double y;
            int i=7;
            int gxmin;
            int gxmax;
            int gymin;
            int gymax;

            if (gametype!=3)
            {
                for (int row=0;row<5;row++)
                {
                    for (int h=-row;h<row+1;h=h+2)
                    {
                        balls[i]._shape.setFillColor(sf::Color(255,0,0));
                        balls[i]._x=pink_x-0.1-2*ball_radius-row*sqrt(3.0)*(ball_radius+DOUBLE_EPSILON);
                        balls[i]._y=pink_y+h*(ball_radius+DOUBLE_EPSILON);
                        balls[i]._order=i+1;
                        i+=1;
                    }
                }
            }
            else
            {
                for (int j=0;j<2;j++)
                {
                    balls[i]._order=i+1;
                    balls[i]._shape.setFillColor(sf::Color(255,0,0));
                    balls[i]._y=black_y;
                    balls[i]._x=rail_thickness+cush_thickness+(black_x-rail_thickness-cush_thickness)*(j+1)/3.;
                    i+=1;
                }
                for (int j=0;j<5;j++)
                {
                    balls[i]._order=i+1;
                    balls[i]._shape.setFillColor(sf::Color(255,0,0));
                    balls[i]._y=black_y;
                    balls[i]._x=black_x+(pink_x-black_x)*(j+1)/6.;
                    i+=1;
                }
                for (int j=0;j<5;j++)
                {
                    balls[i]._order=i+1;
                    balls[i]._shape.setFillColor(sf::Color(255,0,0));
                    balls[i]._y=black_y;
                    balls[i]._x=pink_x+(blue_x-pink_x)*(j+1)/6.;
                    i+=1;
                }
                for (int j=0;j<3;j++)
                {
                    balls[i]._order=i+1;
                    balls[i]._shape.setFillColor(sf::Color(255,0,0));
                    balls[i]._y=black_y;
                    balls[i]._x=blue_x+(blue_x-pink_x)*(j+1)/6.;
                    i+=1;
                }
            }

            for (int i=0;i<22;i++)
            {
                gxmin=int(floor((balls[i]._x-balls[i]._r)/(2.*balls[i]._r)));
                gxmax=int(floor((balls[i]._x+balls[i]._r)/(2.*balls[i]._r)));
                gymin=int(floor((balls[i]._y-balls[i]._r)/(2.*balls[i]._r)));
                gymax=int(floor((balls[i]._y+balls[i]._r)/(2.*balls[i]._r)));

                balls[i]._gpos[0][0]=gxmin;
                balls[i]._gpos[0][1]=gymin;
                balls[i]._gpos[1][0]=gxmax;
                balls[i]._gpos[1][1]=gymin;
                balls[i]._gpos[2][0]=gxmax;
                balls[i]._gpos[2][1]=gymax;
                balls[i]._gpos[3][0]=gxmin;
                balls[i]._gpos[3][1]=gymax;
            }

            for (int i=0;i<22;i++)
            {
                server.serverballs[i]._x=balls[i]._x;
                server.serverballs[i]._y=balls[i]._y;
                server.serverballs[i]._order=balls[i]._order;

                gxmin=int(floor((server.serverballs[i]._x-server.serverballs[i]._r)/(2.*server.serverballs[i]._r)));
                gxmax=int(floor((server.serverballs[i]._x+server.serverballs[i]._r)/(2.*server.serverballs[i]._r)));
                gymin=int(floor((server.serverballs[i]._y-server.serverballs[i]._r)/(2.*server.serverballs[i]._r)));
                gymax=int(floor((server.serverballs[i]._y+server.serverballs[i]._r)/(2.*server.serverballs[i]._r)));

                server.serverballs[i]._gpos[0][0]=gxmin;
                server.serverballs[i]._gpos[0][1]=gymin;
                server.serverballs[i]._gpos[1][0]=gxmax;
                server.serverballs[i]._gpos[1][1]=gymin;
                server.serverballs[i]._gpos[2][0]=gxmax;
                server.serverballs[i]._gpos[2][1]=gymax;
                server.serverballs[i]._gpos[3][0]=gxmin;
                server.serverballs[i]._gpos[3][1]=gymax;
            }

            baulkline.setPrimitiveType(sf::PrimitiveType::LineStrip);
            baulkline.resize(2);
            baulkcircle.setPrimitiveType(sf::PrimitiveType::LineStrip);
            baulkcircle.resize(100);
            baulkline[0].position=sf::Vector2f(_sfac*brown_x,_sfac*(rail_thickness+cush_thickness));
            baulkline[1].position=sf::Vector2f(_sfac*brown_x,(_sfac*raw_height/(1.+panel_ratio))-_sfac*(rail_thickness+cush_thickness));

            for (int i=0;i<100;i++)
            {
                x=brown_x+11.687*sin(pi*i/99);
                y=brown_y+11.687*cos(pi*i/99);
                baulkcircle[i].position=sf::Vector2f(_sfac*x,_sfac*y);
            }

            _shapes.push_back(&baulkline);
            _shapes.push_back(&baulkcircle);

            for (int i=0;i<6;i++)
            {
                _shapes.push_back(&spots[i]);
            }

            cuetraj.setPrimitiveType(sf::PrimitiveType::LineStrip);
            cuetraj2.setPrimitiveType(sf::PrimitiveType::LineStrip);
            obtraj.setPrimitiveType(sf::PrimitiveType::LineStrip);
            ghostball.setRadius(ball_radius*_sfac);
            ghostball.setOrigin(ball_radius*_sfac,ball_radius*_sfac);
            ghostball.setFillColor(sf::Color(255,255,255,100));
            ghostball.setPosition(sf::Vector2f(-100.,-100.));
            ghostball2.setRadius(ball_radius*_sfac);
            ghostball2.setOrigin(ball_radius*_sfac,ball_radius*_sfac);
            ghostball2.setFillColor(sf::Color(255,255,255,100));
            ghostball2.setPosition(sf::Vector2f(-100.,-100.));

            _shapes.push_back(&cuetraj);
            _shapes.push_back(&cuetraj2);
            _shapes.push_back(&obtraj);
            _shapes.push_back(&ghostball);
            _shapes.push_back(&ghostball2);

            for (int i=0;i<22;i++)
            {
                _shapes.push_back(&balls[i]._shape);
            }

            _shapes.push_back(&redarrow._head);
            _shapes.push_back(&redarrow._tail);
            _shapes.push_back(&cue._sprite);

            _shapes.push_back(&pauserect);
            _shapes.push_back(&pausetext);

            _shapes.push_back(&gameoverrect);
            _shapes.push_back(&gameovertext);
            _shapes.push_back(&stats_title);
            for (int i=0;i<stats_text.size();i++)
            {
                _shapes.push_back(&stats_text[i]);
            }
        }
        void update(double dt,sf::Vector2i mouse_pos);
        void scores_update();
};

void GameScreen::scores_update()
{
    //multiplayer.
    textp1score.setString(std::to_string(p1score));
    textp2score.setString(std::to_string(p2score));

    sf::FloatRect bounds;
    bounds=textp1score.getLocalBounds();
    textp1score.setOrigin(sf::Vector2f(bounds.left+0.5*bounds.width,bounds.top+0.5*bounds.height));
    bounds=textp2score.getLocalBounds();
    textp2score.setOrigin(sf::Vector2f(bounds.left+0.5*bounds.width,bounds.top+0.5*bounds.height));

    textp1frames.setString(std::to_string(p1frames));
    textp2frames.setString(std::to_string(p2frames));
    bounds=textp1frames.getLocalBounds();
    textp1frames.setOrigin(sf::Vector2f(bounds.left+0.5*bounds.width,bounds.top+0.5*bounds.height));
    bounds=textp2frames.getLocalBounds();
    textp2frames.setOrigin(sf::Vector2f(bounds.left+0.5*bounds.width,bounds.top+0.5*bounds.height));

    if (isyourturn)
    {
        p1pointer.setFillColor(sf::Color(0,0,0));
        p2pointer.setFillColor(sf::Color(0,0,0,0));
    }
    else
    {
        p1pointer.setFillColor(sf::Color(0,0,0,0));
        p2pointer.setFillColor(sf::Color(0,0,0));
    }

    //gameover stats.
    stats_text[4].setString(std::to_string(p1frames));
    stats_text[5].setString(std::to_string(p2frames));
    stats_text[7].setString(std::to_string(p1_highbreak));
    stats_text[8].setString(std::to_string(p2_highbreak));
    stats_text[10].setString(std::to_string(p1_centuries));
    stats_text[11].setString(std::to_string(p2_centuries));

    if (p1frames>p2frames) {stats_text[4].setFont(_boldfont);stats_text[5].setFont(_thinfont);}
    else if (p1frames<p2frames) {stats_text[4].setFont(_thinfont);stats_text[5].setFont(_boldfont);}
    else {stats_text[4].setFont(_thinfont);stats_text[5].setFont(_thinfont);}

    if (p1_highbreak>p2_highbreak) {stats_text[7].setFont(_boldfont);stats_text[8].setFont(_thinfont);}
    else if (p1_highbreak<p2_highbreak) {stats_text[7].setFont(_thinfont);stats_text[8].setFont(_boldfont);}
    else {stats_text[7].setFont(_thinfont);stats_text[8].setFont(_thinfont);}

    if (p1_centuries>p2_centuries) {stats_text[10].setFont(_boldfont);stats_text[11].setFont(_thinfont);}
    else if (p1_centuries<p2_centuries) {stats_text[10].setFont(_thinfont);stats_text[11].setFont(_boldfont);}
    else {stats_text[10].setFont(_thinfont);stats_text[11].setFont(_thinfont);}

    for (int i=0;i<stats_text.size();i++)
    {
        bounds=stats_text[i].getLocalBounds();
        stats_text[i].setOrigin(sf::Vector2f(bounds.left+0.5*bounds.width,bounds.top+0.5*bounds.height));
    }
}

void GameScreen::update(double dt,sf::Vector2i mouse_pos)
{
    //activate or deactive text and buttons.

    if (gameover || !(fabs(pauserect.getPosition().y+_sfac*raw_height)<0.001))
    {
        _inputboxes[0]._isactive=false;
        _buttons[2]._isactive=false;
        _buttons[3]._isactive=false;
    }
    else
    {
        _inputboxes[0]._isactive=true;
        _buttons[2]._isactive=true;
        _buttons[3]._isactive=true;
    }

    if (!gameover)
    {
        //concession packets.
        if (_buttons[2]._target=="Concedeframe1")
        {
            _buttons[2]._target="Concedeframe";
            if (gametype<2)
            {
                packet.clear();
                packetId=4;
                packet << packetId;
                socket.send(packet);
            }
        }
        else if (_buttons[3]._target=="Concedematch1")
        {
            _buttons[3]._target="Concedematch";
            if (gametype<2)
            {
                packet.clear();
                packetId=5;
                packet << packetId;
                socket.send(packet);
            }
        }

        typing=false;
        for (int i=0;i<_inputboxes.size();i++)
        {
            if (_inputboxes[i]._shape.getOutlineColor()==_inputboxes[i]._outlinecolour2)
            {
                typing=true;
                _inputboxes[i]._t+=dt;
                if (int(_inputboxes[i]._cursor.getFillColor().a)==255 && _inputboxes[i]._t>_inputboxes[i]._ton)
                {
                    _inputboxes[i]._t=0.;
                    sf::Color c=_inputboxes[i]._cursor.getFillColor();
                    _inputboxes[i]._cursor.setFillColor(sf::Color(c.r,c.g,c.b,0));
                }
                else if (int(_inputboxes[i]._cursor.getFillColor().a)==0 && _inputboxes[i]._t>_inputboxes[i]._toff)
                {
                    _inputboxes[i]._t=0.;
                    sf::Color c=_inputboxes[i]._cursor.getFillColor();
                    _inputboxes[i]._cursor.setFillColor(sf::Color(c.r,c.g,c.b,255));
                }
            }
        }

        if (typing)
        {
            if (sf::Keyboard::isKeyPressed(sf::Keyboard::Enter))
            {
                //send the message if not empty.
                if (_inputboxes[0]._input.length()>0)
                {
                    //non-empty message. Send packet as string.
                    packet.clear();
                    packetId=3;
                    packet << packetId << _inputboxes[0]._input;
                    socket.send(packet);

                    _inputboxes[0]._input="";
                    _inputboxes[0]._t=0.;
                    _inputboxes[0]._cursorpos=0;
                    sf::Color c=_inputboxes[0]._cursor.getFillColor();
                    _inputboxes[0]._cursor.setFillColor(sf::Color(c.r,c.g,c.b,255));
                    sf::Vector2f rectpos=_inputboxes[0]._cursor.getPosition();
                    _inputboxes[0]._text.setString(_inputboxes[0]._input.substr(0,_inputboxes[0]._cursorpos));
                    sf::FloatRect bounds=_inputboxes[0]._text.getGlobalBounds();
                    _inputboxes[0]._cursor.setPosition(sf::Vector2f(bounds.left+bounds.width+5.*_inputboxes[0]._abscursorthickness,rectpos.y));
                    _inputboxes[0]._text.setString(_inputboxes[0]._input);
                }
            }
        }

        change=false;
        //listen for packets.
        if (gametype<2)
        {
            packet.clear();
            if (socket.receive(packet)==sf::Socket::Done)
            {
                packet >> packetId;

                if (packetId==0)
                {
                    //display cue trajectory prediction.
                    packet >> power >> cue._angle >> cue._offset >> cue._theta >> cue._alpha >> balls[0]._x >> balls[0]._y;
                    change=true;
                }
                else if (packetId==1)
                {
                    packet >> resultsize;
                    result.clear();
                    for (int i=0;i<resultsize;i++)
                    {
                        for (int j=0;j<66;j++)
                        {
                            packet >> temp[j];
                        }
                        result.push_back(temp);
                    }
                    t=0;
                    done=false;
                }
                else if (packetId==2)
                {
                    //whose turn it is.
                    packet >> isyourturn >> isfoul >> ismiss >> placing_white >> isredon >> isfreeball >> gameover;
                    packet >> p1score >> p2score >> p1frames >> p2frames >> p1_highbreak >> p2_highbreak >> p1_centuries >> p2_centuries;
                    if (done) {scores_update();}
                    if (gameover) {scores_update();}
                }
                else if (packetId==3)
                {
                    //received log msg from server.
                    std::string msgname;
                    std::string msg;
                    packet >> msgname >> msg;

                    std::rotate(logstrings.begin(),logstrings.begin()+1,logstrings.end());

                    for (int i=0;i<numlines;i++)
                    {
                        logtext[i].setString(logstrings[i]);
                    }

                    sf::FloatRect bounds;
                    float x=_inputboxes[0]._shape.getSize().x;
                    std::string total=msgname+": "+msg;
                    int c=0;
                    while (total.length()!=0)
                    {
                        c+=1;
                        logtext[numlines-1].setString(total.substr(0,c));
                        logstrings[numlines-1]=total.substr(0,c);
                        bounds=logtext[numlines-1].getLocalBounds();
                        if (bounds.width>x)
                        {
                            if ((total.substr(c-2,1)).compare(" "))
                            {
                                //hyphen.
                                logtext[numlines-1].setString(total.substr(0,c-2)+"-");
                                logstrings[numlines-1]=total.substr(0,c-2)+"-";
                                total.erase(0,c-2);
                                c=0;
                            }
                            else
                            {
                                logtext[numlines-1].setString(total.substr(0,c-1));
                                logstrings[numlines-1]=total.substr(0,c-1);
                                total.erase(0,c-1);
                                c=0;
                            }
                            //recycle the line.
                            std::rotate(logstrings.begin(),logstrings.begin()+1,logstrings.end());

                            for (int i=0;i<numlines;i++)
                            {
                                logtext[i].setString(logstrings[i]);
                            }
                        }
                        if (c==total.length())
                        {
                            break;
                        }
                    }
                }
            }
        }

        if (int(floor(t*100./framerate))<result.size())
        {
            deltat=t*100./framerate-int(floor(t*100./framerate));
            for (int i=0;i<22;i++)
            {
                if (int(floor(t*100./framerate))==int(result.size()-1))
                {
                    //dont interpolate.
                    balls[i]._x=result[result.size()-1][i*3];
                    balls[i]._y=result[result.size()-1][i*3+1];
                    balls[i]._z=result[result.size()-1][i*3+2];
                }
                else
                {
                    balls[i]._x=result[int(floor(t*100./framerate))][i*3]+deltat*(result[int(floor(t*100./framerate))+1][i*3]-result[int(floor(t*100./framerate))][i*3]);
                    balls[i]._y=result[int(floor(t*100./framerate))][i*3+1]+deltat*(result[int(floor(t*100./framerate))+1][i*3+1]-result[int(floor(t*100./framerate))][i*3+1]);
                    balls[i]._z=result[int(floor(t*100./framerate))][i*3+2]+deltat*(result[int(floor(t*100./framerate))+1][i*3+2]-result[int(floor(t*100./framerate))][i*3+2]);
                }
            }
        }
        if (!done && result.size()>0)
        {
            t=t+1;
            if (int(floor(t*100./framerate))>=result.size())
            {
                done=true;
                for (int i=0;i<22;i++)
                {
                    balls[i]._x=result[result.size()-1][i*3];
                    balls[i]._y=result[result.size()-1][i*3+1];
                    balls[i]._z=result[result.size()-1][i*3+2];
                }

                if (isyourturn) {change=true;}

                if (gametype<2)
                {
                    scores_update();
                }
                else if (gametype==2)
                {
                    //singleplayer vs. AI.
                }
                else if (gametype==3)
                {
                    change=true;
                    //solo lineup.
                    //check if correct ball hit.
                    if (server.ball_hit_order.size()==0) {gameover=true; return;}
                    if (server.ball_potted_order.size()==0) {gameover=true; return;}

                    if (isredon)
                    {
                        for (int i=0;i<server.ball_hit_order.size();i++)
                        {
                            if (server.ball_hit_order[i]<8) {gameover=true; return;}
                        }
                        //check pots.
                        for (int i=0;i<server.ball_potted_order.size();i++)
                        {
                            if (server.ball_potted_order[i]<8) {gameover=true; return;}
                        }
                        //add scores.
                        p1score+=server.ball_potted_order.size();
                        p1_highbreak=p1score;
                    }
                    else
                    {
                        if (server.ball_hit_order.size()>1) {gameover=true; return;}
                        if (server.ball_hit_order[0]!=nom_colour_order) {gameover=true; return;}
                        if (server.ball_potted_order.size()>1) {gameover=true; return;}
                        if (server.ball_potted_order[0]!=nom_colour_order) {gameover=true; return;}
                        p1score+=server.ball_potted_order[0];
                        p1_highbreak=p1score;
                    }
                    isredon=!isredon;

                    stats_text[4].setString("N/A");
                    stats_text[5].setString("N/A");
                    stats_text[7].setString(std::to_string(p1_highbreak));
                    stats_text[8].setString("N/A");
                    if (p1_highbreak>=100) {stats_text[10].setString("1");}
                    stats_text[11].setString("N/A");

                    sf::FloatRect bounds;
                    for (int i=0;i<stats_text.size();i++)
                    {
                        bounds=stats_text[i].getLocalBounds();
                        stats_text[i].setOrigin(sf::Vector2f(bounds.left+0.5*bounds.width,bounds.top+0.5*bounds.height));
                    }
                    server.respot();
                    for (int i=1;i<7;i++)
                    {
                        balls[i]._x=server.serverballs[i]._x;
                        balls[i]._y=server.serverballs[i]._y;
                        balls[i]._z=ball_radius;
                    }

                    //scoreboard stuff.
                    textp1score.setString(std::to_string(p1score));
//                    textp2score.setString(std::to_string(p2score));

                    bounds=textp1score.getLocalBounds();
                    textp1score.setOrigin(sf::Vector2f(bounds.left+0.5*bounds.width,bounds.top+0.5*bounds.height));
//                    bounds=textp2score.getLocalBounds();
//                    textp2score.setOrigin(sf::Vector2f(bounds.left+0.5*bounds.width,bounds.top+0.5*bounds.height));

//                    textp1frames.setString(std::to_string(p1frames));
//                    textp2frames.setString(std::to_string(p2frames));
//                    bounds=textp1frames.getLocalBounds();
//                    textp1frames.setOrigin(sf::Vector2f(bounds.left+0.5*bounds.width,bounds.top+0.5*bounds.height));
//                    bounds=textp2frames.getLocalBounds();
//                    textp2frames.setOrigin(sf::Vector2f(bounds.left+0.5*bounds.width,bounds.top+0.5*bounds.height));
                }
            }
        }

        if (ispaused)
        {
            pos=pauserect.getPosition();
            ds=sqrt(fabs(pos.y)*2.*pausea)*dt;
            pauserect.setPosition(sf::Vector2f(0.,fmin(0.,pos.y+ds)));
            pos=_buttons[0]._shape.getPosition();
            _buttons[0]._shape.setPosition(sf::Vector2f(pos.x,fmin(0.5*_sfac*raw_height,pos.y+ds)));
            _buttons[0]._text.setPosition(sf::Vector2f(pos.x,fmin(0.5*_sfac*raw_height,pos.y+ds)));
            pos=pausetext.getPosition();
            pausetext.setPosition(sf::Vector2f(pos.x,fmin(0.15*_sfac*raw_height,pos.y+ds)));
        }
        else
        {
            pos=pauserect.getPosition();
            ds=-sqrt(fabs(pos.y+_sfac*raw_height)*2.*pausea)*dt;
            pauserect.setPosition(sf::Vector2f(0.,fmax(-_sfac*raw_height,pos.y+ds)));
            pos=_buttons[0]._shape.getPosition();
            _buttons[0]._shape.setPosition(sf::Vector2f(pos.x,fmax(-0.5*_sfac*raw_height,pos.y+ds)));
            _buttons[0]._text.setPosition(sf::Vector2f(pos.x,fmax(-0.5*_sfac*raw_height,pos.y+ds)));
            pos=pausetext.getPosition();
            pausetext.setPosition(sf::Vector2f(pos.x,fmax(-0.85*_sfac*raw_height,pos.y+ds)));
        }

        if (placing_white)
        {
            touching=false;
            for (int i=1;i<22;i++)
            {
                if (balls[0]._potted==false)
                {
                    if (sqrt(pow(balls[0]._x-balls[i]._x,2)+pow(balls[0]._y-balls[i]._y,2))<2*ball_radius)
                    {
                        balls[0]._shape.setFillColor(sf::Color(200,0,0));
                        touching=true;
                        break;
                    }
                }
            }
            if (touching==false)
            {
                balls[0]._shape.setFillColor(sf::Color(255,255,255));
            }
        }

        //get inputs.
        if (sf::Keyboard::isKeyPressed(user_controls["Pause game"]))
        {
            if (!ispausepressed)
            {
                ispausepressed=true;
                ispaused=!ispaused;
            }
        }
        else
        {
            ispausepressed=false;
        }

        if (!typing && isyourturn && done && fabs(pauserect.getPosition().y+_sfac*raw_height)<0.001)
        {
            //get user input when it is their turn.
            if (placing_white)
            {
                if (sf::Keyboard::isKeyPressed(user_controls["Move ball up"]))
                {
                    //move ball up.
                    change=true;
                    if (gametype!=3)
                    {
                        if (sqrt(pow(balls[0]._x-brown_x,2.)+pow(balls[0]._y+0.1*(dt/0.01)-brown_y,2.))<11.687)
                        {
                            balls[0]._y+=0.1*(dt/0.01);
                        }
                    }
                    else
                    {
                        balls[0]._y=std::min(balls[0]._y+0.1*(dt/0.01),rail_thickness+cush_thickness+table_width-ball_radius);
                    }
                }
                if (sf::Keyboard::isKeyPressed(user_controls["Move ball down"]))
                {
                    //move ball down.
                    change=true;
                    if (gametype!=3)
                    {
                        if (sqrt(pow(balls[0]._x-brown_x,2.)+pow(balls[0]._y-0.1*(dt/0.01)-brown_y,2.))<11.687)
                        {
                            balls[0]._y-=0.1*(dt/0.01);
                        }
                    }
                    else
                    {
                        balls[0]._y=std::max(balls[0]._y-0.1*(dt/0.01),rail_thickness+cush_thickness+ball_radius);
                    }
                }
                if (sf::Keyboard::isKeyPressed(user_controls["Move ball left"]))
                {
                    //move ball left.
                    change=true;
                    if (gametype!=3)
                    {
                        balls[0]._x=std::max(balls[0]._x-0.1*(dt/0.01),brown_x);
                    }
                    else
                    {
                        balls[0]._x=std::max(balls[0]._x-0.1*(dt/0.01),rail_thickness+cush_thickness+ball_radius);
                    }
                }
                if (sf::Keyboard::isKeyPressed(user_controls["Move ball right"]))
                {
                    //move ball right.
                    change=true;
                    if (gametype!=3)
                    {
                        if (sqrt(pow(balls[0]._x+0.1*(dt/0.01)-brown_x,2.)+pow(balls[0]._y-brown_y,2.))<11.687)
                        {
                            balls[0]._x+=0.1*(dt/0.01);
                        }
                    }
                    else
                    {
                        balls[0]._x=std::min(balls[0]._x+0.1*(dt/0.01),rail_thickness+cush_thickness+table_length-ball_radius);
                    }
                }
                if (sf::Keyboard::isKeyPressed(user_controls["Place ball"]))
                {
                    //place ball.
                    change=true;
                    if (!touching)
                    {
                        placing_white=false;
                    }
                    if (gametype==3)
                    {
                        server.serverballs[0]._x=balls[0]._x;
                        server.serverballs[0]._y=balls[0]._y;
                    }
                }
            }
            else
            {
                if (sf::Keyboard::isKeyPressed(user_controls["Aim left"]))
                {
                    //move aim left.
                    change=true;
                    cue._angle-=(pi/180.)*(dt/0.01);
                }
                if (sf::Keyboard::isKeyPressed(user_controls["Aim right"]))
                {
                    //move aim right.
                    change=true;
                    cue._angle+=(pi/180.)*(dt/0.01);
                }
                if (sf::Keyboard::isKeyPressed(user_controls["Precise aim left"]))
                {
                    //precise aim left.
                    change=true;
                    cue._angle-=(0.015*pi/180.)*(dt/0.01);
                }
                if (sf::Keyboard::isKeyPressed(user_controls["Precise aim right"]))
                {
                    //precise aim right.
                    change=true;
                    cue._angle+=(0.015*pi/180.)*(dt/0.01);
                }
                if (sf::Keyboard::isKeyPressed(user_controls["Increase power"]))
                {
                    //increase power.
                    change=true;
                    power=fmin(power+0.5*(dt/0.01),100.5);
                }
                if (sf::Keyboard::isKeyPressed(user_controls["Decrease power"]))
                {
                    //decrease power.
                    change=true;
                    power=fmax(power-0.5*(dt/0.01),0.);
                }
                if (sf::Keyboard::isKeyPressed(user_controls["Increase cue elevation"]))
                {
                    //increase cue elevation.
                    change=true;
                    cue._alpha=fmin(cue._alpha+(0.1*pi/180.)*(dt/0.01),0.5*pi);
                }
                if (sf::Keyboard::isKeyPressed(user_controls["Decrease cue elevation"]))
                {
                    //decrease cue elevation.
                    change=true;
                    cue._alpha=fmax(cue._alpha-(0.1*pi/180.)*(dt/0.01),0.);
                }
                if (sf::Keyboard::isKeyPressed(user_controls["Offset cue tip up"]))
                {
                    //offset cue tip up.
                    change=true;
                    dist=sqrt(pow(cue._offset*sin(cue._theta),2.)+pow(cue._offset*cos(cue._theta)+dt*ball_radius,2.));
                    if (dist<ball_radius-0.005/in2m)
                    {
                        cue._theta=atan2(cue._offset*sin(cue._theta),cue._offset*cos(cue._theta)+dt*ball_radius);
                        cue._offset=dist;
                    }
                }
                if (sf::Keyboard::isKeyPressed(user_controls["Offset cue tip down"]))
                {
                    //offset cue tip down.
                    change=true;
                    dist=sqrt(pow(cue._offset*sin(cue._theta),2.)+pow(cue._offset*cos(cue._theta)-dt*ball_radius,2.));
                    if (dist<ball_radius-0.005/in2m)
                    {
                        cue._theta=atan2(cue._offset*sin(cue._theta),cue._offset*cos(cue._theta)-dt*ball_radius);
                        cue._offset=dist;
                    }
                }
                if (sf::Keyboard::isKeyPressed(user_controls["Offset cue tip left"]))
                {
                    //offset cue tip left.
                    change=true;
                    dist=sqrt(pow(cue._offset*sin(cue._theta)-dt*ball_radius,2.)+pow(cue._offset*cos(cue._theta),2.));
                    if (dist<ball_radius-0.005/in2m)
                    {
                        cue._theta=atan2(cue._offset*sin(cue._theta)-dt*ball_radius,cue._offset*cos(cue._theta));
                        cue._offset=dist;
                    }
                }
                if (sf::Keyboard::isKeyPressed(user_controls["Offset cue tip right"]))
                {
                    //offset cue tip right.
                    change=true;
                    dist=sqrt(pow(cue._offset*sin(cue._theta)+dt*ball_radius,2.)+pow(cue._offset*cos(cue._theta),2.));
                    if (dist<ball_radius-0.005/in2m)
                    {
                        cue._theta=atan2(cue._offset*sin(cue._theta)+dt*ball_radius,cue._offset*cos(cue._theta));
                        cue._offset=dist;
                    }
                }
                if (sf::Keyboard::isKeyPressed(user_controls["Strike cueball"]))
                {
                    //strike cueball.
                    cue._speed=1.2*power;
                    cue.perturb();
                    cue.shot();
                    balls[0]._vx=cue._ballv*sin(cue._angle);
                    balls[0]._vy=cue._ballv*cos(cue._angle);
                    balls[0]._xspin=cue._ballparspin*sin(cue._angle)+cue._ballperspin*cos(cue._angle);
                    balls[0]._yspin=cue._ballparspin*cos(cue._angle)-cue._ballperspin*sin(cue._angle);
                    balls[0]._rspin=cue._ballrspin;

                    if (gametype<2)
                    {
                        result.clear();
                        packet.clear();
                        packetId=2;
                        packet << packetId << balls[0]._vx << balls[0]._vy << balls[0]._xspin << balls[0]._yspin << balls[0]._rspin;
                        socket.send(packet);
                        done=false;
                        return;
                    }

                    if (gametype==3)
                    {
                        //solo lineup.
                        server.serverballs[0]._vx=balls[0]._vx;
                        server.serverballs[0]._vy=balls[0]._vy;
                        server.serverballs[0]._xspin=balls[0]._xspin;
                        server.serverballs[0]._yspin=balls[0]._yspin;
                        server.serverballs[0]._rspin=balls[0]._rspin;

                        server.result=server.simulate(server.serverballs,server.servercushions);
                        result=server.result;
                        t=0;
                        done=false;
                        return;
                    }
                }
                if (sf::Mouse::isButtonPressed(sf::Mouse::Left))
                {
                    //ball on.
                    if (!isredon || isfreeball)
                    {
                        for (int i=0;i<6;i++)
                        {
                            dist=sqrt(pow(mouse_pos.x-nomballs[i]._shape.getPosition().x,2.)+pow(mouse_pos.y-nomballs[i]._shape.getPosition().y,2.));
                            if (dist<nomballs[i]._shape.getRadius())
                            {
                                change=true;
                                nom_colour_order=i+2;
                                break;
                            }
                        }
                    }
                }
            }
        }

        for (int i=0;i<22;i++)
        {
            balls[i]._shape.setPosition(sf::Vector2f(_sfac*balls[i]._x,_sfac*raw_height/(1.+panel_ratio)-_sfac*balls[i]._y));
        }
        if (change && !placing_white)
        {
            sf::Color c;
            for (int i=0;i<100;i++)
            {
                c=power_bar[99-i].getFillColor();
                if (i<int(floor(power))) {power_bar[99-i].setFillColor(sf::Color(c.r,c.g,c.b,255));}
                else {power_bar[99-i].setFillColor(sf::Color(c.r,c.g,c.b,0));}
            }

            spin_dot_pos=elevation_ball.getPosition();
            elevation_pointer.setPosition(sf::Vector2f(spin_dot_pos.x+elevation_ball.getRadius()*cos(cue._alpha),spin_dot_pos.y-elevation_ball.getRadius()*sin(cue._alpha)));
            elevation_pointer.setRotation(-cue._alpha*180./pi);

            if (int(cue._alpha*180./pi)<10)
            {
                elevation_display.setString("0"+std::to_string(int(cue._alpha*180./pi))+"°");
            }
            else
            {
                elevation_display.setString(std::to_string(int(cue._alpha*180./pi))+"°");
            }

            spin_selector_pos=spin_selector.getPosition();
            spin_dot.setPosition(sf::Vector2f(spin_selector_pos.x+cue._offset*sin(cue._theta)*spin_selector.getRadius()/ball_radius,spin_selector_pos.y-cue._offset*cos(cue._theta)*spin_selector.getRadius()/ball_radius));

            cue._sprite.setPosition(sf::Vector2f((balls[0]._x-2.*sin(cue._angle))*dfactor,window_height-dfactor*(balls[0]._y-2.*cos(cue._angle))));
            cue._sprite.setRotation((cue._angle+0.5*pi)*180./pi);

            redarrow._head.setPosition(sf::Vector2f(-1000.,-1000.));
            redarrow._tail.setPosition(sf::Vector2f(-1000.,-1000.));

            cue._speed=1.2*power;
            cue.shot();
            balls[0]._vx=cue._ballv*sin(cue._angle);
            balls[0]._vy=cue._ballv*cos(cue._angle);
            balls[0]._xspin=cue._ballparspin*sin(cue._angle)+cue._ballperspin*cos(cue._angle);
            balls[0]._yspin=cue._ballparspin*cos(cue._angle)-cue._ballperspin*sin(cue._angle);
            balls[0]._rspin=cue._ballrspin;

            predict=trajectory(balls,cushions);

            cuetraj.clear();
            cuetraj2.clear();
            obtraj.clear();

            cuetraj.resize(predict[0].size()/2);
            cuetraj2.resize(1+predict[1].size()/2);
            obtraj.resize(predict[2].size()/2);

            if (predict[0].size()>1)
            {
                ghostball.setPosition(sf::Vector2f(dfactor*predict[0][predict[0].size()-2],window_height-dfactor*predict[0][predict[0].size()-1]));
                cuetraj2[0].position=sf::Vector2f(dfactor*predict[0][predict[0].size()-2],window_height-dfactor*predict[0][predict[0].size()-1]);
            }
            else
            {
                ghostball.setPosition(sf::Vector2f(-100.,-100.));
            }

            if (predict[1].size()>1)
            {
                ghostball2.setPosition(sf::Vector2f(dfactor*predict[1][predict[1].size()-2],window_height-dfactor*predict[1][predict[1].size()-1]));
            }
            else
            {
                ghostball2.setPosition(sf::Vector2f(-100.,-100.));
            }

            for (int i=0;i<predict[0].size()/2;i++)
            {
                cuetraj[i].position=sf::Vector2f(dfactor*predict[0][2*i],window_height-dfactor*predict[0][2*i+1]);
            }
            for (int i=1;i<1+predict[1].size()/2;i++)
            {
                cuetraj2[i].position=sf::Vector2f(dfactor*predict[1][2*i-2],window_height-dfactor*predict[1][2*i-1]);
            }
            for (int i=0;i<predict[2].size()/2;i++)
            {
                obtraj[i].position=sf::Vector2f(dfactor*predict[2][2*i],window_height-dfactor*predict[2][2*i+1]);
            }

            if (ghostball2.getPosition().x>=0.)
            {
                //contacts a ball directly.
                spin_dot_pos=ghostball.getPosition();
                dist=999.;
                int contacted=0;
                for (int i=1;i<22;i++)
                {
                    if (sqrt(pow(spin_dot_pos.x-_sfac*balls[i]._x,2.)+pow(spin_dot_pos.y-_sfac*raw_height/(1.+panel_ratio)+_sfac*balls[i]._y,2.))<dist)
                    {
                        dist=sqrt(pow(spin_dot_pos.x-_sfac*balls[i]._x,2.)+pow(spin_dot_pos.y-_sfac*raw_height/(1.+panel_ratio)+_sfac*balls[i]._y,2.));
                        contacted=balls[i]._order;
                    }
                }

                if (contacted<=7)
                {
                    //hits a colour.
                    nom_colour_order=contacted;
                }
            }

            if (isredon)
            {
                //colour choice.
                c=nomtext.getColor();
                nomtext.setColor(sf::Color(c.r,c.g,c.b,100));
                for (int i=0;i<6;i++)
                {
                    nomballs[i]._shape.setFillColor(sf::Color(nomballs[i]._colour.r,nomballs[i]._colour.g,nomballs[i]._colour.b,100));
                    nomballs[i]._shape.setOutlineThickness(0.);
                }
                ballontitle.setString("Ball on - Red");
            }
            else
            {
                c=nomtext.getColor();
                nomtext.setColor(sf::Color(c.r,c.g,c.b,255));
                for (int i=0;i<6;i++)
                {
                    nomballs[i]._shape.setFillColor(sf::Color(nomballs[i]._colour.r,nomballs[i]._colour.g,nomballs[i]._colour.b,255));
                    nomballs[i]._shape.setOutlineThickness(0.);
                }
                if (nom_colour_order>=2 && nom_colour_order<=7)
                {
                    nomballs[nom_colour_order-2]._shape.setOutlineThickness(nomballs[nom_colour_order-2]._abslinethickness);
                }
                if (isfreeball)
                {
                    ballontitle.setString("Ball on - Free ball");
                }
                else
                {
                    ballontitle.setString("Ball on - Colour");
                }
            }

            if (gametype==0 || gametype==1)
            {
                //send the info packet.
                packetId=1;
                packet.clear();

                packet << packetId << power << cue._angle << cue._offset << cue._theta << cue._alpha << balls[0]._x << balls[0]._y;
                socket.send(packet);
            }
        }
        if (placing_white || !done)
        {
            if (placing_white && done)
            {
                redarrow._head.setPosition(sf::Vector2f(dfactor*balls[0]._x,window_height-dfactor*(balls[0]._y+balls[0]._r*2.)));
                redarrow._tail.setPosition(sf::Vector2f(dfactor*balls[0]._x,window_height-dfactor*(balls[0]._y+balls[0]._r*3.)));
            }
            cue._sprite.setPosition(sf::Vector2f(-1000.,-1000.));
            ghostball.setPosition(sf::Vector2f(-1000.,-1000.));
            ghostball2.setPosition(sf::Vector2f(-1000.,-1000.));
            cuetraj.clear();
            cuetraj2.clear();
            obtraj.clear();
        }
    }
    else
    {
        //game over!
        if (fabs(pauserect.getPosition().y+_sfac*raw_height)<0.001)
        {
            pos=gameoverrect.getPosition();
            ds=sqrt(fabs(pos.y)*2.*pausea)*dt;
            gameoverrect.setPosition(sf::Vector2f(0.,fmin(0.,pos.y+ds)));
            pos=_buttons[1]._shape.getPosition();
            _buttons[1]._shape.setPosition(sf::Vector2f(pos.x,fmin(0.75*_sfac*raw_height,pos.y+ds)));
            _buttons[1]._text.setPosition(sf::Vector2f(pos.x,fmin(0.75*_sfac*raw_height,pos.y+ds)));
            pos=gameovertext.getPosition();
            gameovertext.setPosition(sf::Vector2f(pos.x,fmin(0.15*_sfac*raw_height,pos.y+ds)));

            for (int i=0;i<stats_text.size();i++)
            {
                pos=stats_text[i].getPosition();
                stats_text[i].setPosition(sf::Vector2f(pos.x,fmin(stats_text_y[i]+_sfac*raw_height,pos.y+ds)));
            }
        }
        else
        {
            pos=pauserect.getPosition();
            ds=-sqrt(fabs(pos.y+_sfac*raw_height)*2.*pausea)*dt;
            pauserect.setPosition(sf::Vector2f(0.,fmax(-_sfac*raw_height,pos.y+ds)));
            pos=_buttons[0]._shape.getPosition();
            _buttons[0]._shape.setPosition(sf::Vector2f(pos.x,fmax(-0.5*_sfac*raw_height,pos.y+ds)));
            _buttons[0]._text.setPosition(sf::Vector2f(pos.x,fmax(-0.5*_sfac*raw_height,pos.y+ds)));
            pos=pausetext.getPosition();
            pausetext.setPosition(sf::Vector2f(pos.x,fmax(-0.85*_sfac*raw_height,pos.y+ds)));
        }
    }
}


#endif // PLAY-GAME-SCREEN_H_INCLUDED
