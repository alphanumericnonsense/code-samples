# foobar with google

I received this challenge one day while asking too many question on google or something.  I completed everything.  The questions were mostly combinatorial optimization.  Google didn't reach out, so maybe I suck, but it was fun to learn a few things.  I'm scraping this from a notebook I kept, but I didn't record much, so here are some of the statements I found online that matched up with my solutions, probably for the harder ones where I had to do some work.  There's a good list of questions [here](https://github.com/tssovi/google-foobar-challenge/tree/master).

The problems problems below touch on Markov chains, path-finding, maxflow/mincut, and Polya enumeration.

- `markov.py`.  **Doomsday Fuel**.  Making fuel for the LAMBCHOP’s reactor core is a tricky process because of the exotic matter involved. It starts as raw ore, then during processing, begins randomly changing between forms, eventually reaching a stable form. There may be multiple stable forms that a sample could ultimately reach, not all of which are useful as fuel.

  Commander Lambda has tasked you to help the scientists increase fuel creation efficiency by predicting the end state of a given ore sample. You have carefully studied the different structures that the ore can take and which transitions it undergoes. It appears that, while random, the probability of each structure transforming is fixed. That is, each time the ore is in 1 state, it has the same probabilities of entering the next state (which might be the same state). You have recorded the observed transitions in a matrix. The others in the lab have hypothesized more exotic forms that the ore can become, but you haven’t seen all of them.

  Write a function solution(m) that takes an array of array of non-negative integers representing how many times that state has gone to the next state and return an array of integers for each terminal state giving the exact probabilities of each terminal state, represented as the numerator for each state, then the denominator for all of them at the end and in simplest form. The matrix is at most 10 by 10. It is guaranteed that no matter which state the ore is in, there is a path from that state to a terminal state. That is, the processing will always eventually end in a stable state. The ore starts in state 0. The denominator will fit within a signed 32-bit integer during the calculation, as long as the fraction is simplified regularly.

  For example, consider the matrix m:
  ```
  [
    [0,1,0,0,0,1],  # s0, the initial state,goes to s1 and s5 with equal probability
    [4,0,0,3,2,0],  # s1 can become s0, s3, or s4, but with different probabilities
    [0,0,0,0,0,0],  # s2 is terminal, and unreachable (never observed in practice)
    [0,0,0,0,0,0],  # s3 is terminal
    [0,0,0,0,0,0],  # s4 is terminal
    [0,0,0,0,0,0],  # s5 is terminal
  ]
  ```
  So, we can consider different paths to terminal states, such as:
  ```
  s0 -> s1 -> s3
  s0 -> s1 -> s0 -> s1 -> s0 -> s1 -> s4
  s0 -> s1 -> s0 -> s5
  ```
  Tracing the probabilities of each, we find that:
  ```
  s2 has probability 0
  s3 has probability 3/14
  s4 has probability 1/7
  s5 has probability 9/14
  ```
  So, putting that together, and making a common denominator, gives an answer in the form of [s2.numerator, s3.numerator, s4.numerator, s5.numerator, denominator] which is [0, 3, 2, 9, 14].

- `dijkstra.py`.  **Prepare the Bunnies' Escape**.  You’re awfully close to destroying the LAMBCHOP doomsday device  and freeing Commander Lambda’s bunny prisoners, but once they’re free of the prison blocks, the bunnies are going to need to escape Lambda’s  space station via the escape pods as quickly as possible. Unfortunately, the halls of the space station are a maze of corridors and dead ends  that will be a deathtrap for the escaping bunnies. Fortunately,  Commander Lambda has put you in charge of a remodeling project that will give you the opportunity to make things a little easier for the  bunnies. Unfortunately (again), you can’t just remove all obstacles  between the bunnies and the escape pods - at most you can remove one  wall per escape pod path, both to maintain structural integrity of the  station and to avoid arousing Commander Lambda’s suspicions.

   You have maps of parts of the space station, each starting at a  prison exit and ending at the door to an escape pod. The map is  represented as a matrix of 0s and 1s, where 0s are passable space and 1s are impassable walls. The door out of the prison is at the top left  (0,0) and the door into an escape pod is at the bottom right (w-1,h-1).

   Write a function `solution(map)` that generates the length of the shortest path from the prison door to  the escape pod, where you are allowed to remove one wall as part of your remodeling plans. The path length is the total number of nodes you pass through, counting both the entrance and exit nodes. The starting and  ending positions are always passable (0). The map will always be  solvable, though you may or may not need to remove a wall. The height  and width of the map can be from 2 to 20. Moves can only be made in  cardinal directions; no diagonal moves are allowed.

- `floyd_warshall.py`.  **Running With Bunnies**You and your rescued bunny prisoners need to get out of this collapsing  death trap of a space station - and fast! Unfortunately, some of the  bunnies have been weakened by their long imprisonment and can't run very fast. Their friends are trying to help them, but this escape would go a lot faster if you also pitched in. The defensive bulkhead doors have  begun to close, and if you don't make it through in time, you'll be  trapped! You need to grab as many bunnies as you can and get through the bulkheads before they close.

   The time it takes to move from your starting point to all  of the bunnies and to the bulkhead will be given to you in a square  matrix of integers. Each row will tell you the time it takes to get to  the start, first bunny, second bunny, ..., last bunny, and the bulkhead  in that order. The order of the rows follows the same pattern (start,  each bunny, bulkhead). The bunnies can jump into your arms, so picking  them up is instantaneous, and arriving at the bulkhead at the same time  as it seals still allows for a successful, if dramatic, escape. (Don't  worry, any bunnies you don't pick up will be able to escape with you  since they no longer have to carry the ones you did pick up.) You can  revisit different spots if you wish, and moving to the bulkhead doesn't  mean you have to immediately leave - you can move to and from the  bulkhead to pick up additional bunnies if time permits.

   In addition to spending time traveling between bunnies,  some paths interact with the space station's security checkpoints and  add time back to the clock. Adding time to the clock will delay the  closing of the bulkhead doors, and if the time goes back up to 0 or a  positive number after the doors have already closed, it triggers the  bulkhead to reopen. Therefore, it might be possible to walk in a circle  and keep gaining time: that is, each time a path is traversed, the same  amount of time is used or added.

   Write a function of the form solution(times, time_limit)  to calculate the most bunnies you can pick up and which bunnies they  are, while still escaping through the bulkhead before the doors close  for good. If there are multiple sets of bunnies of the same size, return the set of bunnies with the lowest prisoner IDs (as indexes) in sorted  order. The bunnies are represented as a sorted list by prisoner ID, with the first bunny being 0. There are at most 5 bunnies, and time_limit is a non-negative integer that is at most 999.

   For instance, in the case of

   ```
   [
       [0, 2, 2, 2, -1],  # 0 = Start
       [9, 0, 2, 2, -1],  # 1 = Bunny 0
       [9, 3, 0, 2, -1],  # 2 = Bunny 1
       [9, 3, 2, 0, -1],  # 3 = Bunny 2
       [9, 3, 2, 2,  0],  # 4 = Bulkhead
   ]
   ```

   and a time limit of 1, the five inner  array rows designate the starting point, bunny 0, bunny 1, bunny 2, and  the bulkhead door exit respectively. You could take the path:

   Start End Delta Time Status:

   ```
   -   0     -    1 Bulkhead initially open
   0   4    -1    2
   4   2     2    0
   2   4    -1    1
   4   3     2   -1 Bulkhead closes
   3   4    -1    0 Bulkhead reopens; you and the bunnies exit
    
   ```
   
   With this solution, you would pick up bunnies 1 and 2.  This is the best combination for this space station hallway, so the  answer is [1, 2].
   
- `edmonds_karp.py`.  **Escape Pods**.  You've blown up the LAMBCHOP doomsday device and broken the bunnies out of  Lambda's prison - and now you need to escape from the space station as  quickly and as orderly as possible! The bunnies have all gathered in  various locations throughout the station, and need to make their way  towards the seemingly endless amount of escape pods positioned in other  parts of the station. You need to get the numerous bunnies through the  various rooms to the escape pods. Unfortunately, the corridors between  the rooms can only fit so many bunnies at a time. What's more, many of  the corridors were resized to accommodate the LAMBCHOP, so they vary in  how many bunnies can move through them at a time.

   Given the starting room numbers of the groups of bunnies,  the room numbers of the escape pods, and how many bunnies can fit  through at a time in each direction of every corridor in between, figure out how many bunnies can safely make it to the escape pods at a time at peak.

   Write a function solution(entrances, exits, path) that  takes an array of integers denoting where the groups of gathered bunnies are, an array of integers denoting where the escape pods are located,  and an array of an array of integers of the corridors, returning the  total number of bunnies that can get through at each time step as an  int. The entrances and exits are disjoint and thus will never overlap.  The path element path [A] [B] = C describes that the corridor going from A to B can fit C bunnies at each time step.  There are at most 50 rooms  connected by the corridors and at most 2000000 bunnies that will fit at a time.

   For example, if you have:

   ```
   entrances = [0, 1]
   exits = [4, 5]
   path =
   [
       [0, 0, 4, 6, 0, 0],  # Room 0: Bunnies
       [0, 0, 5, 2, 0, 0],  # Room 1: Bunnies
       [0, 0, 0, 0, 4, 4],  # Room 2: Intermediate room
       [0, 0, 0, 0, 6, 6],  # Room 3: Intermediate room
       [0, 0, 0, 0, 0, 0],  # Room 4: Escape pods
       [0, 0, 0, 0, 0, 0],  # Room 5: Escape pods
   ]
   ```

   Then in each time step, the following might happen:

   ```
   0 sends 4/4 bunnies to 2 and 6/6 bunnies to 3
   1 sends 4/5 bunnies to 2 and 2/2 bunnies to 3
   2 sends 4/4 bunnies to 4 and 4/4 bunnies to 5
   3 sends 4/6 bunnies to 4 and 4/6 bunnies to 5
   ```

   So, in total, 16 bunnies could make it to the escape  pods at 4 and 5 at each time step.  (Note that in this example, room 3  could have sent any variation of 8 bunnies to 4 and 5, such as 2/6 and  6/6, but the final solution remains the same.)

- `polya.py`.  **Disorderly Escape** Oh no! You've managed to free the bunny prisoners and escape Commander  Lambdas exploding space station, but her team of elite starfighters has  flanked your ship. If you dont jump to hyperspace, and fast, youll be  shot out of the sky! Problem is, to avoid detection by galactic law enforcement, Commander  Lambda planted her space station in the middle of a quasar quantum flux  field. In order to make the jump to hyperspace, you need to know the  configuration of celestial bodies in the quadrant you plan to jump  through. In order to do *that*, you need to figure out how many  configurations each quadrant could possibly have, so that you can pick  the optimal quadrant through which youll make your jump.

   There's something important to note about quasar quantum  flux fields' configurations: when drawn on a star grid, configurations  are considered equivalent by grouping rather than by order. That is, for a given set of configurations, if you exchange the position of any two  columns or any two rows some number of times, youll find that all of  those configurations are equivalent in that way - in grouping, rather  than order.

   Write a function solution(w, h, s) that takes 3 integers  and returns the number of unique, non-equivalent configurations that can be found on a star grid w blocks wide and h blocks tall where each  celestial body has s possible states.

   Equivalency is defined as above: any two star grids with  each celestial body in the same state where the actual order of the rows and columns do not matter (and can thus be freely swapped around).

   Star grid standardization means that the width and height of the grid will always be between 1 and 12, inclusive.

   And while there are a variety of celestial bodies in each  grid, the number of states of those bodies is between 2 and 20,  inclusive.

   The answer can be over 20 digits long, so return it as a  decimal string.  The intermediate values can also be large, so you will  likely need to use at least 64-bit integers.

   For example, consider w=2, h=2, s=2. We have a 2x2 grid where each celestial body is either in state 0 (for instance, silent) or state 1 (for instance, noisy). We can examine which grids are equivalent by swapping rows and columns.

   ```
   00
   00
   ```

   In the above configuration, all celestial  bodies are "silent" - that is, they have a state of 0 - so any swap of  row or column would keep it in the same state.

   ```
   00 00 01 10
   01 10 00 00
   ```
   
   1 celestial body is emitting noise - that  is, has a state of 1 - so swapping rows and columns can put it in any of the 4 positions.  All four of the above configurations are equivalent.

   ```
   00 11
   11 00    
   ```
   
   2 celestial bodies are emitting noise  side-by-side.  Swapping columns leaves them unchanged, and swapping rows simply moves them between the top and bottom.  In both, the *groupings* are the same: one row with two bodies in state 0, one row with two bodies in state 1, and two columns with one of each state.
   
   ```
   01 10
   01 10
   ```
   
   2 noisy celestial bodies adjacent  vertically. This is symmetric to the side-by-side case, but it is  different because there's no way to transpose the grid.
   
   ```
   01 10
   10 01    
   ```
   
   2 noisy celestial bodies diagonally.  Both have 2 rows and 2 columns that have one of each state, so they are  equivalent to each other.
   
   ```
   01 10 11 11
   11 11 01 10
   ```
   
   3 noisy celestial bodies, similar to the case where only one of four is noisy.
   
   ```
   11
   11    
   ```
   
   4 noisy celestial bodies.
   
   There are 7 distinct, non-equivalent grids in total, so solution(2, 2, 2) would return 7.
   
- `secret_message.py`.  Vigenère cipher with gmail handle as a fun ending to successfully completing all the tasks.
