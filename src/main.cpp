#include <iostream>
#include <fstream>
#include <thread>
#include "objects.h"

int main()
{
    std::cout.precision(std::numeric_limits<double>::max_digits10);
    sf::ContextSettings settings;
    settings.antialiasingLevel=8;
    sf::RenderWindow window(sf::VideoMode(sf::VideoMode::getDesktopMode().width,sf::VideoMode::getDesktopMode().height),"Snooker Game",sf::Style::Fullscreen,settings);

    //shaders.
    sf::Shader shader;
    if (!shader.loadFromFile(_vertexShaderFile,_fragmentShaderFile))
    {
        std::cout << "Error. Shader not loaded." << std::endl;
    }

    std::string custom_input;
    auto it=user_controls.begin();
    std::ifstream file(_userConfigFile);
    if (file.is_open())
    {
        while (getline(file,custom_input))
        {
            it->second=StringToKey(custom_input);
            it++;
        } file.close();
    }

    std::ifstream cuefile(_userCueConfigFile);
    if (cuefile.is_open())
    {
        while (getline(cuefile,custom_input))
        {
            cuetexturefile=custom_input;
        } cuefile.close();
    }

    double sdiff=sf::VideoMode::getDesktopMode().width-dfactor*raw_width;

    sf::View view;
    view.reset(sf::FloatRect(-0.5*sdiff,0.,float(sf::VideoMode::getDesktopMode().width),float(sf::VideoMode::getDesktopMode().height)));
    window.setView(view);

    //prepare matrix.
    //KEEP THIS HERE!!! ESSENTIAL FOR TRAJECTORY PREDICTION!!!
    for (int i=0;i<44;i++)
    {
        M_(i,i)=1/ball_mass;
    }
    M_(44,44)=0;
    M_(45,45)=0;

    std::vector<GameState*> states;
    states.push_back(new TitleScreen());
    sf::Vector2f rectsize;
    sf::Vector2f rectpos;
    sf::FloatRect bounds;
    sf::Event event;
    sf::Vector2i mouse_pos;
    bool hastyped=false;
    int typingindex=0;

    std::string name;
    std::string ip;
    unsigned short port;
    std::string localip;
    unsigned short localport;

    Server server;
    localip=server.serverIp.toString();
    localport=server.port;

    std::thread server_thread(&Server::executionThread,std::ref(server));
    server_thread.detach();

    const int bufferSize=25;
    std::vector<double> timeBuffer(bufferSize,1.0/framerate);
    int bufferHead=0;
    double averageTime=1.0/framerate;

    sf::Time elapsed;
    sf::Clock waitingClock;
    sf::Clock clock;

    while (window.isOpen())
    {
        if (states.size()==0) {break;}

        if (!states.back()->_isWaitingForInput)
        {
            states.back()->handleButtonPress(window);

            //handle redirects.
            if (states.back()->_shouldChangeState==true)
            {
                //reset for next time.
                states.back()->_shouldChangeState=false;

                if (states.back()->_stateTarget=="TitleScreen") {states.push_back(new TitleScreen());}
                else if (states.back()->_stateTarget=="Singleplayer") {states.push_back(new SingleplayerScreen());}
                else if (states.back()->_stateTarget=="SingleplayerAI") {states.push_back(new GameScreen(dfactor,2,"",50000,name));}
                else if (states.back()->_stateTarget=="SingleplayerLineup") {states.push_back(new GameScreen(dfactor,3,"",50000,name));}
                else if (states.back()->_stateTarget=="Multiplayer")
                {
                    states.push_back(new MultiplayerScreen(dfactor,localip,std::to_string(localport)));

                    if (server._connectionError)
                    {
                        //display permanent error message for off-line.
                        MultiplayerScreen* temp=dynamic_cast<MultiplayerScreen*>(static_cast<GameState*>(states.back()));
                        temp->displayError("No internet connection!",false);
                    }
                }
                else if (states.back()->_stateTarget=="MultiplayerHost")
                {
                    if (!server._connectionError)
                    {
                        server.resetServer();
                        MultiplayerScreen* temp=dynamic_cast<MultiplayerScreen*>(static_cast<GameState*>(states.back()));
                        server.framesbestof=temp->_framesBestOf;
                        try
                        {
                            states.push_back(new GameScreen(dfactor,0,localip,localport,std::string(states.back()->_inputboxes[2]._text.getString())));
                        }
                        catch (...) {}
                    }
                }
                else if (states.back()->_stateTarget=="MultiplayerJoin")
                {
                    try
                    {
                        ip=states.back()->_inputboxes[0]._text.getString();
                        port=std::stoi(std::string(states.back()->_inputboxes[1]._text.getString()));
                        states.push_back(new GameScreen(dfactor,1,ip,port,std::string(states.back()->_inputboxes[2]._text.getString())));
                    }
                    catch (...)
                    {
                        if (!server._connectionError)
                        {
                            //display an error message.
                            MultiplayerScreen* temp=dynamic_cast<MultiplayerScreen*>(static_cast<GameState*>(states.back()));
                            temp->displayError("Could not connect to host!");
                        }
                    }
                }
                else if (states.back()->_stateTarget=="Options") {states.push_back(new OptionsScreen());}
                else if (states.back()->_stateTarget=="Controls") {states.push_back(new ControlScreen());}
                else if (states.back()->_stateTarget=="Changecue") {states.push_back(new ChangeCueScreen());}
                else if (states.back()->_stateTarget=="Quit") {delete states.back(); states.pop_back();}

                if (states.size()>0) {states.back()->_shouldUpdate=true;}

                continue;
            }

            states.back()->handleInputClick(window);
            typingindex=states.back()->_typingIndex;
        }

        hastyped=false;
        //deal with user input.
        while (window.pollEvent(event))
        {
            switch (event.type)
            {
                case sf::Event::Closed:
                    window.close();
                    break;
                case sf::Event::TextEntered:
                    if (!states.back()->_isWaitingForInput && states.back()->_isTyping && (event.text.unicode<127 && event.text.unicode>31))
                    {
                        hastyped=true;
                        states.back()->_inputboxes[typingindex]._input.insert(size_t(states.back()->_inputboxes[typingindex]._cursorpos),std::string(sf::String(event.text.unicode)));
                        states.back()->_inputboxes[typingindex]._cursorpos+=1;
                    }
                case sf::Event::KeyPressed:
                    if (!states.back()->_isWaitingForInput && states.back()->_isTyping && event.key.code==sf::Keyboard::Backspace)
                    {
                        hastyped=true;
                        if (states.back()->_inputboxes[typingindex]._cursorpos!=0)
                        {
                            states.back()->_inputboxes[typingindex]._input.erase(states.back()->_inputboxes[typingindex]._cursorpos-1,1);
                            states.back()->_inputboxes[typingindex]._cursorpos-=1;
                        }
                    }
                    else if (!states.back()->_isWaitingForInput && states.back()->_isTyping && event.key.code==sf::Keyboard::Left)
                    {
                        hastyped=true;
                        states.back()->_inputboxes[typingindex]._cursorpos=std::max(0,states.back()->_inputboxes[typingindex]._cursorpos-1);
                    }
                    else if (!states.back()->_isWaitingForInput && states.back()->_isTyping && event.key.code==sf::Keyboard::Right)
                    {
                        hastyped=true;
                        states.back()->_inputboxes[typingindex]._cursorpos=std::min(int(states.back()->_inputboxes[typingindex]._input.length()),states.back()->_inputboxes[typingindex]._cursorpos+1);
                    }
                    else if (!states.back()->_isWaitingForInput && states.back()->_isTyping && event.key.code==sf::Keyboard::Delete)
                    {
                        hastyped=true;
                        if (states.back()->_inputboxes[typingindex]._cursorpos!=states.back()->_inputboxes[typingindex]._input.length())
                        {
                            states.back()->_inputboxes[typingindex]._input.erase(states.back()->_inputboxes[typingindex]._cursorpos,1);
                        }
                    }
                    else if (states.back()->_isWaitingForInput && states.back()->_controlindex!=-1)
                    {
                        int controlindex=states.back()->_controlindex;
                        user_controls[states.back()->_buttons[controlindex]._target]=event.key.code;
                        states.back()->_buttons[controlindex]._text.setString(KeyToString(event.key.code));
                        bounds=states.back()->_buttons[controlindex]._text.getLocalBounds();
                        states.back()->_buttons[controlindex]._text.setOrigin(sf::Vector2f(int(bounds.left+0.5*bounds.width),int(bounds.top+0.5*bounds.height)));

                        //update the config file with new controls.
                        std::ofstream file(_userConfigFile,std::ofstream::out | std::ofstream::trunc);
                        for (auto thing:user_controls) {file << KeyToString(thing.second) << "\n";}
                        file.close();

                        //reset the button colour.
                        states.back()->_buttons[controlindex]._shape.setFillColor(states.back()->_buttons[controlindex]._colour1);
                        states.back()->_buttons[controlindex]._shape.setOutlineColor(states.back()->_buttons[controlindex]._outlinecolour1);
                        states.back()->_buttons[controlindex]._text.setFillColor(states.back()->_buttons[controlindex]._textcolour1);

                        states.back()->_shouldUpdate=true;

                        states.back()->_isWaitingForInput=false;
                        states.back()->_controlindex=-1;
                        break;
                    }
            }
        }

        if (!states.back()->_isWaitingForInput && hastyped)
        {
            states.back()->_inputboxes[typingindex]._t=0.;
            sf::Color c=states.back()->_inputboxes[typingindex]._cursor.getFillColor();
            states.back()->_inputboxes[typingindex]._cursor.setFillColor(sf::Color(c.r,c.g,c.b,255));
            rectpos=states.back()->_inputboxes[typingindex]._cursor.getPosition();
            states.back()->_inputboxes[typingindex]._text.setString(states.back()->_inputboxes[typingindex]._input.substr(0,states.back()->_inputboxes[typingindex]._cursorpos));
            bounds=states.back()->_inputboxes[typingindex]._text.getGlobalBounds();
            states.back()->_inputboxes[typingindex]._cursor.setPosition(sf::Vector2f(bounds.left+bounds.width+5.*states.back()->_inputboxes[typingindex]._abscursorthickness,rectpos.y));
            states.back()->_inputboxes[typingindex]._text.setString(states.back()->_inputboxes[typingindex]._input);

            bounds=states.back()->_inputboxes[typingindex]._cursor.getGlobalBounds();
            rectpos=states.back()->_inputboxes[typingindex]._shape.getPosition();
            rectsize=states.back()->_inputboxes[typingindex]._shape.getSize();
            if (bounds.left+bounds.width>rectpos.x+0.5*rectsize.x)
            {
                states.back()->_inputboxes[typingindex]._input.resize(states.back()->_inputboxes[typingindex]._input.length()-1);
                if (states.back()->_inputboxes[typingindex]._cursorpos>states.back()->_inputboxes[typingindex]._input.length())
                {
                    states.back()->_inputboxes[typingindex]._cursorpos=states.back()->_inputboxes[typingindex]._input.length();
                }
                rectpos=states.back()->_inputboxes[typingindex]._cursor.getPosition();
                states.back()->_inputboxes[typingindex]._text.setString(states.back()->_inputboxes[typingindex]._input.substr(0,states.back()->_inputboxes[typingindex]._cursorpos));
                bounds=states.back()->_inputboxes[typingindex]._text.getGlobalBounds();
                states.back()->_inputboxes[typingindex]._cursor.setPosition(sf::Vector2f(bounds.left+bounds.width+5.*states.back()->_inputboxes[typingindex]._abscursorthickness,rectpos.y));
                states.back()->_inputboxes[typingindex]._text.setString(states.back()->_inputboxes[typingindex]._input);
            }
        }

        mouse_pos=sf::Mouse::getPosition(window);
        mouse_pos.x=int(mouse_pos.x-0.5*sdiff);
        if (!states.back()->_isWaitingForInput && states.back()->_shouldUpdate)
        {
            states.back()->handleButtonCallbacks();
            states.back()->update(1.0/framerate,mouse_pos);
        }

        states.back()->updateFpsCounter();
        states.back()->render(window);

        elapsed=clock.restart();
        averageTime+=(elapsed.asSeconds()-timeBuffer[bufferHead])/bufferSize;
        timeBuffer[bufferHead]=elapsed.asSeconds();
        bufferHead=(bufferHead+1)%bufferSize;
        framerate=std::min(int(1.0/averageTime),max_framerate);
        framerate=std::max(framerate,min_framerate);

        waitingClock.restart();
        while (waitingClock.getElapsedTime().asSeconds()<1.0/max_framerate-elapsed.asSeconds()) {};
    }

    return 0;
}
