#ifndef CUSTOM-RECORDER_H_INCLUDED
#define CUSTOM-RECORDER_H_INCLUDED

class CustomRecorder: public sf::SoundRecorder
{
    private:
        unsigned short port;
        sf::IpAddress ipAddress;
        sf::TcpSocket socket; //blocking by default.
        sf::Packet packet;

        virtual bool onStart()
        {
            return true;
        }

        virtual bool onProcessSamples(const sf::Int16* samples, std::size_t sampleCount)
        {
            packet.clear();
            packet.append(samples, sampleCount*sizeof(sf::Int16));
            socket.send(packet);

            return true;
        }

        virtual void onStop() {}
    public:
        void setTarget(sf::IpAddress newIp, unsigned short newPort)
        {
            ipAddress=newIp;
            port=newPort;
            try
            {
                socket.connect(ipAddress,port);
            }
            catch (...) {}
        }
};

#endif // CUSTOM-RECORDER_H_INCLUDED
