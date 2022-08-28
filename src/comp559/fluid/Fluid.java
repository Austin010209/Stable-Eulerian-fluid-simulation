package comp559.fluid;

import java.util.LinkedList;
import java.util.List;

import javax.swing.JPanel;
import javax.swing.border.TitledBorder;
import javax.vecmath.Point2f;
import javax.vecmath.Tuple2f;
import javax.vecmath.Vector2f;

import mintools.parameters.BooleanParameter;
import mintools.parameters.DoubleParameter;
import mintools.parameters.IntParameter;
import mintools.swing.VerticalFlowPanel;

/**
 * Eulerian fluid simulation class. 
 * 
 * This follows closely the solver documented in J. Stam GDC 2003, specifically
 * this uses a non staggered grid, with an extra set of grid cells to help enforce
 * boundary conditions
 * 
 * @author kry
 */
public class Fluid {
    
    final static int DIM = 2;
    
    /** 
     * Velocity, non-staggered, packed (first index is the dimension) 
     * that is, U[0][IX(0,0)] is the x velocity in the (0,0) grid location.
     * */ 
    public float[][] U0;
    
    /** temporary velocity variable */
    private float[][] U1;          
       
    /** temperature (packed)*/
    public float[] temperature0;
    
    /** temperature (packed) temporary variable*/
    private float[] temperature1;
    
    private IntParameter Nval = new IntParameter( "grid size", 16, 4, 256 );
    
    /** Number of grid cells (not counting extra boundary cells */
    public int N = 16;
    
    /** Dimension of each grid cell */
    public float dx = 1;
    
    /** time elapsed in the fluid simulation */
    public double elapsed;

    /**
     * Sources of heat and cold
     */
    public List<Source> sources = new LinkedList<Source>();
    
    /**
     * initialize memory
     */
    public void setup() {
        elapsed = 0;
        N = Nval.getValue();        
        dx = 1.0f / N; // we choose the domain size here to be 1 unit square!
        
        int np2s = (N+2)*(N+2);
        U0 = new float[2][np2s];
        U1 = new float[2][np2s];
        temperature0 = new float[np2s];
        temperature1 = new float[np2s];
    }

    /**
     * Compute the index 
     * @param i 
     * @param j 
     * @return the index in the flattened array
     */
    public int IX( int i, int j ) {
//        return i*(N+2) + j;
    	return j*(N+2) + i;
    }
    
    /**
     * Adjusts values in the boundary cells such that interpolation near the boundary provides 
     * the desired values.  The result depends on the specified flag b
     * if b == 0 then only continuity is guaranteed
     * if b == 1 then the value is fixed so that the interpolated quantity will goes to zero at the x boundaries
     * if b == 2 then the value is fixed so that the interpolated quantity will goes to zero at the y boundaries
     * @param b
     * @param x
     */
    public void setBoundary( int b, float[] x ) {
        int i;
        for ( i=1 ; i<=N; i++ ) {
            x[IX(0 ,i)]  = b==1  ? -x[IX(1,i)] : x[IX(1,i)];
            x[IX(N+1,i)] = b==1 ? -x[IX(N,i)] : x[IX(N,i)];
            x[IX(i,0 )]  = b==2  ? -x[IX(i,1)] : x[IX(i,1)];
            x[IX(i,N+1)] = b==2 ? -x[IX(i,N)] : x[IX(i,N)];            
        }
        x[IX(0 ,0 )] = 0.5f*(x[IX(1,0 )]+x[IX(0 ,1)]);
        x[IX(0 ,N+1)] = 0.5f*(x[IX(1,N+1)]+x[IX(0 ,N )]);
        x[IX(N+1,0 )] = 0.5f*(x[IX(N,0 )]+x[IX(N+1,1)]);
        x[IX(N+1,N+1)] = 0.5f*(x[IX(N,N+1)]+x[IX(N+1,N )]);
    }

    
    /** 
     * Gets the velocity at the given point using interpolation 
     * @param x
     * @param vel
     */
    public void getVelocity( Tuple2f x, Tuple2f vel ) {
        getVelocity( x, U0, vel );
    }
    
    /** 
     * Gets the velocity in the provided velocity field at the given point using interpolation
     * @param x
     * @param U
     * @param vel
     */
    private void getVelocity( Tuple2f x, float[][] U, Tuple2f vel ) {
        vel.x = interpolate( x, U[0] );
        vel.y = interpolate( x, U[1] );
    }
    
    /**
     * Interpolates the given scalar field
     * @param x Location in the grid
     * @param s grid of quantities
     * @return interpolated value
     */
    public float interpolate( Tuple2f x, float[] s ) {
    	
    	// TODO: Objective 1: implement bilinear interpolation (try to make this code fast!)
    	float x0 = x.x / dx - 0.5f;
        float y0 = x.y / dx - 0.5f;

        int xfl = (int) x0;
        int xcei = xfl + 1;

        int yfl = (int) y0;
        int ycei = yfl + 1;

        // Edge Cases
        float LL = 0, LR = 0, RL = 0, RR = 0;
        if (xfl < 0) {
        	LL = -1;
        	LR = -1;
        }

        if (xcei > N + 1) {
        	RL = -1;
        	RR = -1;
        }

        if (yfl < 0) {
        	LL = -1;
        	RL = -1;
        }

        if (ycei > N + 1) {
        	LR = -1;
        	RR = -1;
        }


        if (LL == -1) {
        	LL = 0;
        } else {
        	int idx = IX(xfl, yfl);
        	if(idx < s.length && idx > 0) {
        		LL = s[idx];
        	}
        }

        if (LR == -1) {
        	LR = 0;
        } else {
        	int idx = IX(xfl, ycei);
        	if(idx < s.length && idx > 0) {
        		LR = s[idx];
        	}
        }

        if (RL == -1) {
        	RL = 0;
        } else {
        	int idx = IX(xcei, yfl);
        	if(idx < s.length && idx > 0) {
        		RL = s[idx];
        	}
        }

        if (RR == -1) {
        	RR = 0;
        } else {
        	int idx = IX(xcei, ycei);
        	if(idx < s.length && idx > 0) {
        		RR = s[idx];
        	}
        }
        
        float ydist = y0 - yfl;
        float Ydist = 1 - ydist;
        float xdist = x0 - xfl;
        float f0 = Ydist * LL + ydist * LR;
        float f1 = Ydist * RL + ydist * RR;

        return f0 * (1 - xdist) + f1 * xdist;

    }
        
    /** 
     * Performs a simple Forward Euler particle trace using the current velocity field 
     * @param x0 Current particle location
     * @param h	 Time step
     * @param x1 Final particle location
     */
    public void traceParticle( Point2f x0, float h, Point2f x1 ) {
        traceParticle( x0, U0, h, x1 );        
    }
    
    /** 
     * Performs a simple particle trace using Forward Euler.  Up to the caller
     * to provide a positive or negative time step depending if they want to 
     * trace forward (i.e., filaments) or backwards (advection term).
     * Note that this could also be a higher order integration, or adaptive!
     * x1 = x0 + h * U(x0)
     * @param x0   Starting point
     * @param U    Velocity field
     * @param h    Time step
     * @param x1   Resulting point
     */
    private void traceParticle( Point2f x0, float[][] U, float h, Point2f x1 ) {

    	// TODO: Objective 4: Implement the tracing of a particle position in the velocity field.
    	// Use the getVelocity method, which calls the interpolation method (that you need to write)
    	Tuple2f Ux0 = new Vector2f();
    	getVelocity(x0, U, Ux0 );
    	x1.set(Ux0);
    	x1.scale(h);
    	x1.add(x0);
    }
        
    /**
     * Diffuse the given scalar field by the given amount.  Enforce the specified boundary conditions on the result.
     * @param S1   diffused quantities
     * @param S0   initial quantities
     * @param b    boundary conditions 0 = continuity, 1 = goes to zero on x boundary, 2 = goes to zero on y boundaries
     * @param diff diffusion coefficient
     * @param dt   time step
     */
    private void diffuse( float[] S1, float[] S0, int b, float diff, float dt ) {
        
    	// TODO: Objective 3: Implement diffusion on the given grid of quantities
    	float a = dt*diff*N*N;
    	for(int k = 0; k < iterations.getValue(); k++) {
    		for(int i = 1; i <= N; i++) {
    			for(int j = 1; j <= N; j++) {
    				S1[IX(i,j)] = (S0[IX(i,j)] + a*( S1[IX(i-1,j)] + 
    						S1[IX(i+1,j)] + S1[IX(i,j-1)] + S1[IX(i,j+1)] )) / (1+4*a);
    			}
    		}
    		setBoundary(b, S1);
    	}
    }
    
    /**
     * Advects / transports scalar quantities by the given velocities.
     * @param s1 	Final advected quantities
     * @param s0    Current grid of quantities
     * @param U		Velocity field
     * @param dt 	Time step
     */
    public void transport( float[] s1, float[] s0, float[][] U, float dt ) {
        
    	// TODO: Objective 5: Implement implicit advection of quantities by tracing particles backwards in time in the provided velocity field.

    	int xfl, yfl, xcei, ycei;
        float newx, newy, rx0, ry0, rx1, ry1;
        
        Point2f pt = new Point2f();
        Point2f newpt = new Point2f();

        for (int i = 1; i <= N; i++) {
        	for (int j = 1; j <= N; j++) {
	            pt.x = dx * (i + 0.5f);
	            pt.y = dx * (j + 0.5f);
	            traceParticle(pt, U, -dt, newpt);
	            newx = newpt.x / dx - 0.5f;
	            newy = newpt.y / dx - 0.5f;
	            
	            if (newx < 0.5) newx = 0.5f;
	            if (newx > N + 0.5) newx = N + 0.5f;
	            xfl = (int) newx;
	            xcei = xfl + 1;
	
	            if (newy < 0.5) newy = 0.5f;
	            if (newy > N + 0.5) newy = N + 0.5f;
	            yfl = (int) newy;
	            ycei = yfl + 1;
	
	            rx1 = newx - xfl;
	            rx0 = 1 - rx1;
	            ry1 = newy - yfl;
	            ry0 = 1 - ry1;
	            s1[IX(i, j)] = rx0 * (ry0 * s0[IX(xfl, yfl)] + ry1 * s0[IX(xfl, ycei)])
	                    + rx1 * (ry0 * s0[IX(xcei, yfl)] + ry1 * s0[IX(xcei, ycei)]);
          }
        }
        setBoundary(0, s1);
    }
    
    /**
     * Does the Poisson solve to make sure that velocities U respect incompressible flow
     * @param U
     */
    private void project( float[][] U ) {    

    	// TODO: Objective 6: Implement pressure projection on the provided velocity field
    	
    	//p and div in the paper are the variable we create and will use!
    	int size = (N+2)*(N+2);
    	float[] p = new float[size];
    	float[] div = new float[size];
    	float h = (float) (1.0/N);
    	for(int i = 1; i <= N; i++) {
    		for(int j = 1; j <= N; j++) {
    			div[IX(i,j)] = (float) ( -0.5*h*( U[0][IX(i+1,j)] - U[0][IX(i-1,j)] 
    					+ U[1][IX(i,j+1)] - U[1][IX(i,j-1)] ) );
    			p[IX(i,j)] = 0;
    		}
    	}
    	setBoundary(0, div);
    	setBoundary(0, p);
    	
    	for(int k = 0; k < iterations.getValue(); k++) {
    		for(int i = 1; i <= N; i++) {
    			for(int j = 1; j <= N; j++) {
    				p[IX(i,j)] = ( div[IX(i,j)] + p[IX(i-1,j)] + 
    						p[IX(i+1,j)] + p[IX(i,j-1)] + p[IX(i,j+1)] ) / 4;
    			}
    		}
    		setBoundary(0, p);
    	}
    	
    	for(int i = 1; i <= N; i++) {
    		for(int j = 1; j <= N; j++) {
    			U[0][IX(i,j)] -= 0.5*( p[IX(i+1,j)] - p[IX(i-1,j)] )/h;
    			U[1][IX(i,j)] -= 0.5*( p[IX(i,j+1)] - p[IX(i,j-1)] )/h;
    		}
    	}
    	setBoundary(1, U[0]);
    	setBoundary(2, U[1]);
    }
    
    
    /**
     * Adds a force at a given point in the provided velocity field.
     * @param U
     * @param dt
     * @param x
     * @param f
     */
    private void addForce( float[][] U, float dt, Tuple2f x, Tuple2f f ) {
        addSource( U[0], dt, x, f.x );
        addSource( U[1], dt, x, f.y );
    }
    
    /**
     * Adds some time step scaled amount to the provided scalar field.  
     * Use bilinear interpolation to distribute the amount to the 4 closest cells.
     * @param S			quantity field to modify
     * @param dt		time step
     * @param x			position
     * @param amount	amount
     */
    private void addSource( float[] S, float dt, Tuple2f x, float amount ) {

    	// TODO: Objective 2: add a "source" to the provided quantity field.  
    	// Use bilinear interpolation (similar to your interpolate method) to distribute the amount.
    	// Note that this is used by mouse interaction and temperature forces on the velocity field (through addForce)
    	// as well as for heat sources and sinks (i.e., user created points in the grid).
 
    	double X = x.x / dx - 0.5f;
        double Y = x.y / dx - 0.5f;
        int xfl = (int) X;
        int xcei = xfl + 1;
        int yfl = (int) Y;
        int ycei = yfl + 1;

        boolean xfl_bound = (xfl < 0.5 || xfl >= N + 0.5);
        boolean xcei_bound = (xcei < 0.5 || xcei >= N + 0.5);
        boolean yfl_bound = (yfl < 0.5 || yfl >= N + 0.5);
        boolean ycei_bound = (ycei < 0.5 || ycei >= N + 0.5);

        
        double xdist = xcei - X;
        double ydist = ycei - Y;
        double ydist1 = 1.0 - ydist;
        double xfl_amt = xdist * amount;
        double xceil_amt = amount - xfl_amt;

        
        int indx1 = IX(xfl, yfl);
        int indx2 = IX(xfl, ycei);
        int indx3 = IX(xcei, yfl);
        int indx4 = IX(xcei, ycei);
        if(indx1 > 0 && indx1 < S.length && !(xfl_bound || yfl_bound)) {
            S[indx1] += dt * ydist * xfl_amt;
        }
        if(indx2 > 0 && indx2 < S.length && !(xfl_bound || ycei_bound)) {
            S[indx2] += dt * ydist1 * xfl_amt;
        }
        if(indx3 > 0  && indx3 < S.length && !(xfl_bound || ycei_bound)) {
            S[indx3] += dt * ydist * xceil_amt;
        }
        if(indx4 > 0 && indx4 < S.length && !(xcei_bound || ycei_bound)) {
            S[indx4] += dt * ydist1 * xceil_amt;
        }
    }
    
    /**
     * Gets the average temperature of the continuum.  (use this in computing buoyancy forces)
     * @return  average temperature
     */
    public double getReferenceTemperature() {
    	int count = 0;
        double referenceTemperature = 0;
        for ( int i = 1; i <= N; i++ ) {
            for ( int j = 1; j <= N; j++ ) {
                referenceTemperature += temperature0[IX(i,j)];
                count++;
            }
        }
        referenceTemperature /= count;
        return referenceTemperature;
    }
    
    /**
     * Applies buoyancy forces to the velocity field due to temperature.
	 * Use Foster and Metaxis [1997] Equation 2:	F_{bv} = beta g_v (T_0 - T_k)  
     * @param U
     * @param dt
     */
    private void addTemperatureForce( float[][] U, float dt ) {       	
        double referenceTemperature = getReferenceTemperature();
        float beta = buoyancy.getFloatValue();
        
    	// TODO: Objective 7: change velocities based on the temperature.  Don't forget to set Boundaries after modifying velocities!
        for(int i = 1; i <= N; i++) {
        	for(int j = 1; j <= N; j++) {
        		U[1][IX(i,j)] += beta * dt * (referenceTemperature - temperature0[IX(i,j)]);
        	}
        }
        setBoundary(2, U[1]);

        
        
    }
    
    /** Worker variables for mouse interaction */
    private Point2f XVX = new Point2f();
    private Point2f Xprev = new Point2f();
    
    /** 
     * Add forcing along a mouse drag vector.
     * The mouse dragging positions are set by a previous call to setMouseMotionPos.
     */
    private void addMouseForce( float[][] U, float dt ) {
        Vector2f f = new Vector2f();
        f.sub( XVX, Xprev );
        float d = Xprev.distance(XVX);
        if ( d < 1e-6 ) return;
        f.scale( mouseForce.getFloatValue() );
        // go along the path of the mouse!
        Point2f x = new Point2f();
        int num = (int) (d/dx + 1);
        for ( int i = 0; i <= num; i++ ) {
            x.interpolate(Xprev,XVX, (float)i / num );
            addForce( U, dt, x, f );
        }
        Xprev.set( XVX );
    }
    
    /**
     * Sets the mouse location in the fluid for doing dragging interaction
     * @param x0
     * @param x1
     */
    public void setMouseMotionPos( Point2f x0, Point2f x1 ) {
        Xprev.set( x0 );
        XVX.set( x1 );
    }
    
    /**
     * Performs the velocity step
     * @param dt
     */
    public void velocityStep( float dt ) {
        float visc = viscosity.getFloatValue();
    	float[][] tmp;
    	if ( velocityDiffuse.getValue() ) {
	        diffuse( U1[0], U0[0], 1, visc, dt );
	        diffuse( U1[1], U0[1], 2, visc, dt );     
	        tmp = U1; U1 = U0; U0 = tmp;
    	}
    	if ( velocityProject.getValue() ) {
    		project ( U0 ); 
    	}
    	if ( velocityAdvect.getValue() ) {
	        transport( U1[0], U0[0], U1, dt );
	        transport( U1[1], U0[1], U1, dt );
	        tmp = U1; U1 = U0; U0 = tmp;
	        setBoundary( 1, U0[0] ); 
	        setBoundary( 2, U0[1] );
    	}    
        if ( velocityProject.getValue() ) {
        	project ( U0 );
        }
    }
    
    // controls for enabling and disabling different steps, possibly useful for debugging!
    private BooleanParameter velocityDiffuse = new BooleanParameter( "velocity step diffuse", true );
    private BooleanParameter velocityProject = new BooleanParameter( "velcity step project", true );
    private BooleanParameter velocityAdvect = new BooleanParameter( "velocity step advect", true );
    private BooleanParameter scalarDiffuse = new BooleanParameter( "scalar step diffuse", true );
    private BooleanParameter scalarAdvect = new BooleanParameter( "scalar step advect", true );
    
    /**
     * performs the scalar step
     * @param dt
     */
    public void scalarStep(float dt) {
    	float[] tmpt;        
    	if ( scalarDiffuse.getValue() ) {
	    	float diff= diffusion.getFloatValue();
	        diffuse( temperature1, temperature0, 0, diff, dt );        
	        tmpt = temperature1; temperature1 = temperature0; temperature0 = tmpt;
    	}
        
    	if ( scalarAdvect.getValue() ) {
	        transport( temperature1, temperature0, U0, dt );
	        tmpt = temperature1; temperature1 = temperature0; temperature0 = tmpt;
	        setBoundary(0, temperature0); 
    	}    	
    }
    
    /**
     * Advances the state of the fluid by one time step
     */
    public void step() {
        float dt = timeStepSize.getFloatValue();
        addMouseForce( U0, dt );
        for ( Source s : sources ) {
            addSource( temperature0, dt, s.location, s.amount );
        }
        addTemperatureForce( U0, dt );
        velocityStep(dt);
        scalarStep(dt);        
        elapsed += dt;
    }
        
    private DoubleParameter viscosity = new DoubleParameter( "viscosity", 1e-6, 1e-8, 1 );
    private DoubleParameter diffusion = new DoubleParameter( "diffusion", 1e-6, 1e-8, 1 );
    private DoubleParameter buoyancy = new DoubleParameter( "bouyancy", 0.1, -1, 1 );
    private IntParameter iterations = new IntParameter( "GS iterations", 30, 0, 1000 );    
    private DoubleParameter mouseForce = new DoubleParameter( "mouse force", 1e2, 1, 1e3 );
    public DoubleParameter timeStepSize = new DoubleParameter( "time step size", 0.1, 0.001, 1 );
    
    /**
     * Get the parameters for the fluid 
     * @return a control panel
     */
    public JPanel getControls() {
        VerticalFlowPanel vfp = new VerticalFlowPanel();

        VerticalFlowPanel vfp1 = new VerticalFlowPanel();
        vfp1.setBorder(new TitledBorder("Fluid properties"));
        vfp1.add( viscosity.getSliderControls(true) );
        vfp1.add( diffusion.getSliderControls(true) );
        vfp1.add( buoyancy.getSliderControls(false) );
        vfp1.add( mouseForce.getSliderControls(true) );
        vfp.add( vfp1.getPanel() );

        VerticalFlowPanel vfp2 = new VerticalFlowPanel();
        vfp2.setBorder(new TitledBorder("Fluid solver properties"));
        vfp2.add( Nval.getSliderControls() );
        vfp2.add( timeStepSize.getSliderControls(true ) ); 
        vfp2.add( iterations.getSliderControls() );
        vfp.add( vfp2.getPanel() );
        
        VerticalFlowPanel vfp3 = new VerticalFlowPanel();
        vfp3.setBorder( new TitledBorder("Disable/Enable steps"));        
		vfp3.add(velocityDiffuse.getControls());
		vfp3.add(velocityProject.getControls());
		vfp3.add(velocityAdvect.getControls());
		vfp3.add(scalarDiffuse.getControls());
		vfp3.add(scalarAdvect.getControls());
        vfp.add( vfp3.getPanel() );
        		
        return vfp.getPanel();
    }
}