

![Release Screenshot 2022 12 07 - 18 03 53 77](https://user-images.githubusercontent.com/53403691/206262466-675fb7b2-7bca-4eab-ae4f-c5cfbd1896f2.png)

# Snooker

This is a basic 2D snooker game. The game contains a functional UI as well as custom controls and the ability to choose your cue design.

To help with aiming, the game includes a shot predictor to show the trajectory of the cueball.

# Supported game modes:
- Singleplayer line up: Place the cueball anywhere and make a break from there
- LAN multiplayer: Play against another person and communicate via live text chat
- (Coming soon) Singleplayer vs. AI: Play against the computer

# Screenshots

![Release Screenshot 2022 12 07 - 18 05 08 55](https://user-images.githubusercontent.com/53403691/206262059-cc5cf94a-88da-4e5d-bf37-9f0c83dd8485.png)

![Release Screenshot 2022 12 07 - 18 03 28 13](https://user-images.githubusercontent.com/53403691/206262071-5c2b0178-998c-4235-8175-9374ac3d83d0.png)

![Release Screenshot 2022 12 07 - 18 10 12 36](https://user-images.githubusercontent.com/53403691/206262189-476eee2d-3acd-44da-9c93-fc596b022fdf.png)


# Example Gameplay

https://user-images.githubusercontent.com/53403691/206261928-a69ba1ab-b702-4fbc-a262-8e9eb1063979.mp4

https://user-images.githubusercontent.com/53403691/206261988-412e2b0c-4cea-4531-9dbc-2e1f696a2d2d.mp4

# Physics

The general motion of the balls is solved analytically rather than using a discrete timestep for greater accuracy and efficiency [1].

However the motion of a ball on the rim of a pocket is solved numerically [2].

To solve the quartic equations required for the simulation, a self-coded implementation of Algorithm 954 was used [3].

Collisions between balls were solved simulataneously using a matrix method [4]. While collisions between discrete pairs of balls could be solved using a simpler method, the matrix method was useful for providing realistic break-off shots at the start of a frame.

# References

[1] W. Leckie and M. Greenspan, "An Event-Based Pool Physics Simulator", Lecture Notes in Computer Science, pp. 247-262, 2006. Available: 10.1007/11922155_19.

[2] M. Bacon, "How balls roll off tables", American Journal of Physics, vol. 73, no. 8, pp. 722-724, 2005. Available: 10.1119/1.1947198.

[3] N. Flocke, "Algorithm 954", ACM Transactions on Mathematical Software, vol. 41, no. 4, pp. 1-24, 2015. Available: 10.1145/2699468.

[4] M. balls, A. Zita, A. Zita, J. Alexiou and s. gerbil, "Multiple colliding balls", Physics Stack Exchange, 2022. [Online]. Available: https://physics.stackexchange.com/questions/296767/multiple-colliding-balls. [Accessed: 28- Aug- 2022].
