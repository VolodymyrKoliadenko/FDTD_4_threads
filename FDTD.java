
import static java.lang.Math.*;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
//import java.util.Calendar;

public class FDTD {

    public static ExecutorService service = Executors.newFixedThreadPool(4);// fixed threads

    public static float[][] callFDTD(int nx, int ny, String method) {
        CountDownLatch cdl;

        int i, j;//for loops
        float x; //for coordinate
        final float lambd = 1064e-9f; //final is GOOD !!!!
        final float dx = lambd / 15;         // max dx=lam/10 !!
        final float dy = dx;
        final float period = 2e-6f;
        final float Q = 1.08f;//1.08
        final float n = 1;//refractive index

        final float prodol = 2 * n * period * period / lambd / Q;

        float[][] Ez = new float[nx + 1][ny + 1];
        float[][] Hx = new float[nx - 1 + 1][ny - 1 + 1];
        float[][] Hy = new float[nx][ny];

        final float omega = (float) (2 * PI / lambd);
        final float dt = dx / 2;
        final float tau = 2e-5f * 999;// light decay
        final float Z = 376.7303f;
        final float s = dt / dx;
        final float k3 = (1 - s) / (1 + s);//for MUR
        final float w = 19e-7f;//19 gauss
        final float alpha = (float) (sin(0.0 / 180 * PI));//radians. АККУРАТНЕЙ с целочисленным делением

        final int begin = 10;// must bee small
        System.out.println("Imeem__Sloy= " + (ny - begin) * dy / prodol * 2);
        final float mod = 0.008f * 2;//max modulation of EPSILON =2*deltaN;

        final float ds = dt * dt / dx / dx;// renormalization  dH=dt/dx/Z;
        float[][] e = new float[nx + 1][ny + 1];
        for (i = 1; i < nx + 1; i++) {
            for (j = 1; j < ny + 1; j++) {
                // EPSILON MATRIX full initialization
                //ПЕРЕНОРМИРОВКА МОДУЛЯЦИИ
                e[i][j] = (float) (ds / (n + ((j < begin) ? 0 : (mod / 2) * (1 + signum(-0.1 + cos(2 * PI * (i - nx / 2.0 + 0.5) * dx / period) * sin(2 * PI * (j - begin) * dy / prodol))))));
            }
        }

        float[][] end = new float[2][nx + 1]; // boundary conditions
        float[][] top = new float[2][ny + 1];
        float[][] bottom = new float[2][ny + 1];
        //long tTime  = Calendar.getInstance().getTimeInMillis();

        final int tmax = (int) (ny * 2.2);
        System.out.println("START CICLE");
        for (int t = 1; t <= tmax; t++) {                    //____nachalo cicla
            float tt = Math.min(t * s + 10, ny - 1);
            //gauss // bez zatuhania          // int/2
            switch (method) {
                case "cos":        //напишы тут функцию вместо 4 копипа
                    for (i = 1; i <= nx - 1; i++) {
                        x = (float) (dx * (i - (float) nx / 2 + 0.5));
                        Ez[i][1] = (float) (exp(-pow(x, 2) / w / w - (t - 1) * dt / tau) * cos((x * alpha + (t - 1) * dt) * omega));//перепишешь лямбді массива
                    }
                    break;
                case "sin":
                    for (i = 1; i <= nx - 1; i++) {
                        x = (float) (dx * (i - (float) nx / 2 + 0.5));
                        Ez[i][1] = (float) (exp(-pow(x, 2) / w / w - (t - 1) * dt / tau) * sin((x * alpha + (t - 1) * dt) * omega));
                    }
                    break;
            }

            for (i = 1; i <= nx; i++) {  // boundary conditions
                end[0][i] = Ez[i][ny - 1];
                end[1][i] = Ez[i][ny];
            }
            System.arraycopy(Ez[1], 0, top[0], 0, ny + 1);
            System.arraycopy(Ez[2], 0, top[1], 0, ny + 1);
            System.arraycopy(Ez[nx - 1], 0, bottom[0], 0, ny + 1);
            System.arraycopy(Ez[nx], 0, bottom[1], 0, ny + 1);

            // main Ez  
            cdl = new CountDownLatch(2);
            NewThreadE first = new NewThreadE(Ez, Hx, Hy, e, 2, nx / 2, tt, cdl);
            NewThreadE second = new NewThreadE(Ez, Hx, Hy, e, nx / 2 + 1, nx - 1, tt, cdl);
            service.execute(first);
            service.execute(second);
            try {
                cdl.await();
            } catch (InterruptedException ex) {
            }

            for (i = 1; i <= nx; i++) {    // boundary conditions
                Ez[i][ny] = end[0][i] + k3 * (end[1][i] - Ez[i][ny - 1]);//end
            }
            for (i = 1; i <= ny; i++) {
                Ez[1][i] = top[1][i] + k3 * (top[0][i] - Ez[2][i]);//verh kray
                Ez[nx][i] = bottom[0][i] + k3 * (bottom[1][i] - Ez[nx - 1][i]);
            }
            //Ez=Arrays.stream(Ez0).map(float[]::clone).toArray(float[][]::new);  VERY SLOW !!! 2.7X
            //Ez=Arrays.copyOf(Ez0, Ez0.length);   FAST
            switch (method) {
                case "cos":
                    for (i = 1; i <= nx - 1; i++) {
                        x = (float) (dx * (i - (float) nx / 2 + 0.5));
                        Ez[i][1] = (float) (exp(-pow(x, 2) / w / w - (t - 1) * dt / tau) * cos((x * alpha + t * dt) * omega));
                    }
                    break;
                case "sin":
                    for (i = 1; i <= nx - 1; i++) {
                        x = (float) (dx * (i - (float) nx / 2 + 0.5));
                        Ez[i][1] = (float) (exp(-pow(x, 2) / w / w - (t - 1) * dt / tau) * sin((x * alpha + t * dt) * omega));
                    }
                    break;
            }

            cdl = new CountDownLatch(2);
            NewThreadH firstH = new NewThreadH(Ez, Hx, Hy, 1, nx / 2, tt, cdl);
            NewThreadH secondH = new NewThreadH(Ez, Hx, Hy, nx/2+1, nx-1, tt, cdl);
            service.execute(firstH);
            service.execute(secondH);
            try {
                cdl.await();
            } catch (InterruptedException ex) {
            }

        }

        //System.out.println((float)(Calendar.getInstance().getTimeInMillis() - tTime)/1000 + " sec");
        int pos = method.equals("cos") ? 0 : 1;
        BasicEx.forFurier[pos] = new double[nx];
        int endF = (int) (ny * 0.95);//0.99
        for (i = 1; i <= nx; i++) {
            BasicEx.forFurier[pos][i - 1] = Ez[i][endF];
            for (j = 1; j <= ny; j++) {
                Ez[i][j] = abs(Ez[i][j]);// ABS
            }
        }
        Hx = null;
        Hy = null;
        e = null;
        return Ez;

    }

}
