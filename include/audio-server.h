#ifndef AUDIO-SERVER_H_INCLUDED
#define AUDIO-SERVER_H_INCLUDED

class AudioServer
{
    public:
        bool running=true;

        sf::IpAddress serverIp;
        unsigned short port;
        sf::TcpListener listener;

        std::array<sf::TcpSocket,2> playersAudioRecorder;
        std::array<sf::TcpSocket,4> spectatorsAudioRecorder;

        std::array<sf::TcpSocket,2> playersAudioReceiver;
        std::array<sf::TcpSocket,4> spectatorsAudioReceiver;

        AudioServer();
        void executionThread();
};

AudioServer::AudioServer()
{
    //intialize audio server.
}

void AudioServer::executionThread()
{
    while (running)
    {
        //do stuff.
    }
}

#endif // AUDIO-SERVER_H_INCLUDED
