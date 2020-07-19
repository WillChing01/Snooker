# Snooker

This is a basic 2D snooker game. The physics is simulated using the techniques described in 'An Event-Based Pool Physics Simulator' by Leckie and Greenspan (2006). Collisions between balls are solved analytically rather than using a discrete timestep to increase accuracy.

The physics still has to be tweaked slightly - sidespin is not right at the moment.

To help with aiming, the game includes a basic shot predictor to show the trajectory of the cueball.

# Controls

- Left or Right arrow - aim the cue
- Up or Down arrow - increase/decrease shot power
- "<" or ">" - aim the cue more precisely
- "[" or "]" - increase/decrease cue elevation
- Spacebar - take the shot

If the white is potted, using the arrows to position it in the 'D' and press Spacebar to place it.

Click within the white cueball on the bottom left to select where to strike the ball to impart spin.

# Things to improve:

- Implement sidespin more accurately
- GUI with menus/options
- Multiplayer function
- Basic AI to practise against
