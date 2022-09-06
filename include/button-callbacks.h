#ifndef BUTTON-CALLBACKS_H_INCLUDED
#define BUTTON-CALLBACKS_H_INCLUDED

//change cue screen.

void scrollLeft(GameState* game_state, std::map<std::string,std::string>* payload)
{
    ChangeCueScreen* state=dynamic_cast<ChangeCueScreen*>(game_state);
    state->_page=std::max(state->_page-1,0);
}

void scrollRight(GameState* game_state, std::map<std::string,std::string>* payload)
{
    ChangeCueScreen* state=dynamic_cast<ChangeCueScreen*>(game_state);
    state->_page=std::min(state->_page+1,(state->totalCues-1)/state->numPerPage);
}

void selectCue(GameState* game_state, std::map<std::string,std::string>* payload)
{
    ChangeCueScreen* state=dynamic_cast<ChangeCueScreen*>(game_state);

    state->selectedCueIndex=std::stoi((*payload)["Select"]);
    cuetexturefile=_cueFilePrefix+"cue"+std::to_string(state->selectedCueIndex)+".png";
    std::ofstream cuefile(_userCueConfigFile,std::ofstream::out | std::ofstream::trunc);
    cuefile << cuetexturefile; cuefile.close();
}

//control screen.

void setDefaultControls(GameState* game_state, std::map<std::string,std::string>* payload)
{
    ControlScreen* state=dynamic_cast<ControlScreen*>(game_state);

    user_controls=default_controls;
    sf::FloatRect bounds;
    for (int j=0;j<default_controls.size();j++)
    {
        state->_buttons[j]._text.setString(KeyToString(user_controls[state->_buttons[j]._target]));
        bounds=state->_buttons[j]._text.getLocalBounds();
        state->_buttons[j]._text.setOrigin(sf::Vector2f(int(bounds.left+0.5*bounds.width),int(bounds.top+0.5*bounds.height)));
    }

    //update the config file with new controls.
    std::ofstream file(_userConfigFile,std::ofstream::out | std::ofstream::trunc);
    for (auto thing:user_controls) {file << KeyToString(thing.second) << "\n";}
    file.close();
}

void setWaitingForControl(GameState* game_state, std::map<std::string,std::string>* payload)
{
    ControlScreen* state=dynamic_cast<ControlScreen*>(game_state);

    int i=std::stoi((*payload)["index"]);
    sf::FloatRect bounds;
    state->_controlindex=i;
    state->_buttons[i]._text.setString("?");
    bounds=state->_buttons[i]._text.getLocalBounds();
    state->_buttons[i]._text.setOrigin(sf::Vector2f(int(bounds.left+0.5*bounds.width),int(bounds.top+0.5*bounds.height)));
    state->_isWaitingForInput=true;
}

#endif // BUTTON-CALLBACKS_H_INCLUDED
