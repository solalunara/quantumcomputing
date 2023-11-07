using quantumbits;
using static quantumbits.complex_utils;

qbit a = new( 0, 1 );
qbit b = false;
qbit c = false;

qbit_group system = new( a, b, c );

//pre-transmission
/*
b.LeftRotation();
c = b ^ c;
b = a ^ b;
a.RightRotation();
*/

system = system.LeftShift( 1 );
system = system.Xor( 1, 2 );
system = system.Xor( 0, 1 );
system = system.RightShift( 0 );


Console.WriteLine( system );

//post-transmission
/*
a.S_Shift();
c = b ^ c;
a = c ^ a;
a.S_Shift();
c.T_Shift();
a = c ^ a;
*/
system = system.SShift( 0 );
system = system.Xor( 1, 2 );
system = system.Xor( 2, 0 );
system = system.SShift( 0 );
system = system.TShift( 2 );
system = system.Xor( 2, 0 );

system.TryDisentangle( out qbit?[] qbits );

Console.WriteLine( system );