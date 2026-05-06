#include <QPainter>
#include <stdio.h>
#include <cmath>
#include <QKeyEvent>
#include "inter_app.h"

#include "window.h"

#define DEFAULT_A -10
#define DEFAULT_B 10
#define DEFAULT_N 10
#define DEFAULT_K 0
#define EPS 1e-16

static double f_0(double x)
{
    return 1.0 + 0.0 * x;
}

static double f_1(double x)
{
    return x;
}

static double f_2(double x)
{
	return pow(x, 2);
}

static double f_3(double x)
{
	return pow(x, 3);
}

static double f_4(double x)
{
	return pow(x, 4);
}

static double f_5(double x)
{
	return exp(x);
}

static double f_6(double x)
{
	return 1/(25*x*x + 1);
}



Window::Window(QWidget *parent) : QWidget(parent)
{
    a = DEFAULT_A;
    b = DEFAULT_B;
    n = DEFAULT_N;
    k = DEFAULT_K;
    func_id = k;
    display = 0;
    scale = 0;
    perturbation = 0;
    max_f = 0.0;
    setFocusPolicy(Qt::StrongFocus);
    change_func();
}

QSize Window::minimumSizeHint() const
{
    return QSize(100, 100);
}

QSize Window::sizeHint() const
{
    return QSize(1000, 1000);
}

int Window::parse_command_line(int argc, char *argv[])
{
    if (argc != 5) {
	        qWarning("Wrong amount of arguments\n");
		return -1;
	}
    if (sscanf(argv[1], "%lf", &a) != 1 ||
        sscanf(argv[2], "%lf", &b) != 1 ||
        b - a < 1.e-6) {
        qWarning("Wrong arguments a or b");
        return -2;
    }
    if (sscanf(argv[3], "%d", &n) != 1 || n <= 0) {
        qWarning("Wrong argument n");
        return -2;
    }
    if (sscanf(argv[4], "%d", &k) != 1 || k < 0 || k > 6) {
        qWarning("Wrong argument k");
        return -2;
    }
    func_id = k;
    change_func();
    return 0;
}

/// change current function for drawing
void Window::change_func()
{
    func_id = (func_id + 1) % 7;

    switch (func_id) {
    case 0:
        f_name = "f (x) = 1";
        f = f_0;
        break;
    case 1:
        f_name = "f (x) = x";
        f = f_1;
        break;
    case 2:
        f_name = "f (x) = x^2";
        f = f_2;
        break;
    case 3:
        f_name = "f (x) = x^3";
        f = f_3;
        break;
    case 4:
        f_name = "f (x) = x^4";
        f = f_4;
        break;
    case 5:
        f_name = "f (x) = e^x";
        f = f_5;
        break;
    case 6:
        f_name = "f (x) = 1/(25*x^2 + 1)";
        f = f_6;
        break;
    }
    update();
}

/// render graph
void Window::paintEvent(QPaintEvent * /* event */)
{
    QPainter painter(this);
    double s = pow(2.0, scale);
	double x1 = 0.0;
	double x2 = 0.0;
	double y1 = 0.0;
	double y2 = 0.0;
    double a_scaled = a / s;
    double b_scaled = b / s;
	double x = 0.0;
    double delta_x = (b_scaled - a_scaled) / n;
	double delta_y = 0.0;
	double max_abs = 0.0;
    double max_y = 0.0;
	double min_y = 0.0;
	double value = 0.0;
	double*X = new double[n];
	double*F = new double[n];
	double *A1 = nullptr;
	double *A2 = nullptr;
	double*extra1 = nullptr;
	double*extra2 = nullptr;
	double step = (b_scaled - a_scaled)/ (n - 1);
	bool is_second_method =  false;
	bool is_first_method =  false;
	bool first_p = true;
	int draw_points = 1000;
	if(n <= 50) { draw_points = n * 10;}
	if (draw_points < 100) draw_points = 100;
	delta_x = (b_scaled - a_scaled) / draw_points;

	for(int i = 0; i < n; i++) {
		X[i] = a_scaled + i * step;
		F[i] = f(X[i]);

		// 7e
		if(i == n/2) {
			F[i] += perturbation * 0.1 * max_f;
		}
	}

	if(n >= 4) {
		A1 = new double[4 * (n - 1)];
		extra1 = new double[3 *n];
		if(BuildingAkima(n, X, F, A1, extra1) == 0) {
			is_first_method = true;
		}
		A2 = new double[4 * (n - 1)];
		extra2 = new double[5 *n];
//		if(BuildingSplines(n, X, F, A2, extra2) == 0) {
//			is_second_method = true;
//		}
	}

	for(double i = a_scaled; i <= b_scaled; i += delta_x) {
		switch(display) {
			case 0: // Akima
			value = f(x);
			if(first_p || value < min_y) { min_y = value; }
			if(first_p || value > max_y) { max_y = value; }
			if(is_first_method == true) {
				value = EvaluationAkima(x, a_scaled, b_scaled, n, X, A1);
				if(value < min_y) { min_y = value; }
				if(value > max_y) { max_y = value; }
			}
			break;

			case 1: // splines
			value = f(x);
			if(first_p || value < min_y) { min_y = value; }
			if(first_p || value > max_y) { max_y = value; }
			if(is_second_method == true) {
//				value = EvaluationSplines(x, a_scaled, b_scaled, n, X, A2);
				if(value < min_y) { min_y = value; }
				if(value > max_y) { max_y = value; }
			}
			break;

			case 2: // both
			value = f(x);
			if(first_p || value < min_y) { min_y = value; }
			if(first_p || value > max_y) { max_y = value; }
			if(is_first_method == true) {
				value = EvaluationAkima(x, a_scaled, b_scaled, n, X, A1);
				if(value < min_y) { min_y = value; }
				if(value > max_y) { max_y = value; }
			}
			if(is_second_method == true) {
//				value = EvaluationSplines(x, a_scaled, b_scaled, n, X, A2);
				if(value < min_y) { min_y = value; }
				if(value > max_y) { max_y = value; }
			}
			break;

			case 3:
			if(is_first_method == true) {
				value = fabs(EvaluationAkima(x, a_scaled, b_scaled, n, X, A1) - f(x));
				if(first_p || value < min_y) { min_y = value; }
				if(first_p || value > max_y) { max_y = value; }
				first_p = false;
			}
			if(is_second_method == true) {
//				value = fabs(EvaluationSplines(x, a_scaled, b_scaled, n, X, A2) - f(x));
				if(value < min_y) { min_y = value; }
				if(value > max_y) { max_y = value; }
				first_p = false;
			}
			break;
		}
		first_p = false;
	}


    for (int k = 0; k < 2; k++) {
        x = a_scaled;
	if(k == 1) { x = b_scaled; }      
        if (display != 3) {
            value = f(x);
            if (first_p || value < min_y) { min_y = value;}
            if (first_p || value > max_y) { max_y = value; }
            first_p = false;
        }
        if (is_first_method == true) {
            if (display != 3) {
                value = EvaluationAkima(x, a_scaled, b_scaled, n, X, A1);
		}
            else {
                value = fabs(EvaluationAkima(x, a_scaled, b_scaled, n, X, A1) - f(x));
		}
            if (first_p || value < min_y) {min_y = value; }
            if (first_p || value > max_y) {max_y = value; }
            first_p = false;
        }
        if (is_second_method == true) {
            if (display != 3) {
//                value = EvaluationSplines(x, a_scaled, b_scaled, n, X, A2);
		}
            else {
//                value = fabs(EvaluationSplines(x, a_scaled, b_scaled, n, X, A2) - f(x));
		}
            if (first_p || value < min_y) { min_y = value; }
            if (first_p || value > max_y) { max_y = value; }
            first_p = false;
        }
    }



	max_abs = fabs(min_y);
	if(fabs(max_y) > fabs(min_y)) { max_abs = fabs(max_y); }
	max_f = max_abs;

	delta_y = 0.05 * (max_y - min_y);
	if(delta_y < EPS) { delta_y = 1.0; }
	min_y -= delta_y;
	max_y += delta_y;


	painter.save();
	painter.translate(0, height());
	painter.scale(width() / (b_scaled - a_scaled), -height() / (max_y - min_y));
	painter.translate(-a_scaled, -min_y);
        QPen pen("black");
        pen.setWidth(0);
        painter.setPen(pen);

	// drawing functions
	switch (display) {
	    case 0: 
	        pen.setColor("red");
	        painter.setPen(pen);
	        DrawingFunction(painter, a_scaled, b_scaled, delta_x);        
        if (is_first_method == true) {
            pen.setColor("black");
		pen.setWidth(0);
            painter.setPen(pen);
            DrawingApproximation(painter, a_scaled, b_scaled, delta_x, n, X, A1, 1);
		pen.setStyle(Qt::DashLine);
		pen.setStyle(Qt::SolidLine);
        }
        break;
        
	    case 1: 
	        pen.setColor("red");
	        painter.setPen(pen);
	        DrawingFunction(painter, a_scaled, b_scaled, delta_x);
	        if (is_second_method == true) {
	            pen.setColor("black");
		pen.setWidth(0);
	            painter.setPen(pen);
	            DrawingApproximation(painter, a_scaled, b_scaled, delta_x, n, X, A2, 2);
		pen.setStyle(Qt::DashLine);
		pen.setStyle(Qt::SolidLine);
	        }
        break;
        
	    case 2: 
	        pen.setColor("red");
	        painter.setPen(pen);
	        DrawingFunction(painter, a_scaled, b_scaled, delta_x);
                if (is_first_method == true) {
            pen.setColor("black");
		pen.setWidth(0);
            painter.setPen(pen);
            DrawingApproximation(painter, a_scaled, b_scaled, delta_x, n, X, A1, 1);
		pen.setStyle(Qt::DashLine);
		pen.setStyle(Qt::SolidLine);
	        }
	         if (is_second_method == true) {
	            pen.setColor("black");
		pen.setWidth(0);
	            painter.setPen(pen);
	            DrawingApproximation(painter, a_scaled, b_scaled, delta_x, n, X, A2, 2);
		pen.setStyle(Qt::DashLine);
		pen.setStyle(Qt::SolidLine);
	        }
        break;
        
	    case 3: 
	        pen.setColor("red");
	        painter.setPen(pen);
	        DrawingFunction(painter, a_scaled, b_scaled, delta_x);
                if (is_first_method == true) {
            pen.setColor("black");
		pen.setWidth(0);
            painter.setPen(pen);
            DrawingError(painter, a_scaled, b_scaled, delta_x, n, X, A1, 1);
		pen.setStyle(Qt::DashLine);
		pen.setStyle(Qt::SolidLine);
	        }
	         if (is_second_method == true) {
	            pen.setColor("black");
		pen.setWidth(0);
	            painter.setPen(pen);
	            DrawingError(painter, a_scaled, b_scaled, delta_x, n, X, A2, 2);
		pen.setStyle(Qt::DashLine);
		pen.setStyle(Qt::SolidLine);
	        }
        break;
    }

    
    x1 = a_scaled, y1 = f(x1);
    for (x2 = x1 + delta_x; x2 <= b_scaled; x2 += delta_x) {
        y2 = f(x2);
        painter.drawLine(QPointF(x1, y1), QPointF(x2, y2));
        x1 = x2; y1 = y2;
    }
    
    pen.setColor("red");
    painter.setPen(pen);
    painter.drawLine(a_scaled, 0, b_scaled, 0);
    painter.drawLine(0, min_y, 0, max_y);
    
    painter.restore();
    
    painter.setPen("blue");
    painter.drawText(10, 20, QString("k=%1, %2").arg(k).arg(f_name));
    painter.drawText(10, 40, QString("max|f| = %1").arg(max_abs));
    painter.drawText(10, 60, QString("scale = %1").arg(scale));
    painter.drawText(10, 80, QString("n = %1").arg(n));
    painter.drawText(10, 100, QString("perturbation p = %1").arg(perturbation));
    painter.drawText(10, 120, QString("mode = %1").arg(display));
    delete[] X;
    delete[] F;
    if (A1)  {delete[] A1; }
    if (extra1) { delete[] extra1; }
    if (A2) { delete[] A2; }
    if (extra2) { delete[] extra2; }

}



// drawing one of seven functions
void Window::DrawingFunction(QPainter &painter, double a, double b, double dx)
{
	QPen pen = painter.pen(); 
    painter.setPen(pen);
    double x1 = a;
	double y1= f(x1);
    for (double x2 = x1 + dx; x2 <= b + dx*0.5; x2 += dx) {
        painter.drawLine(QPointF(x1, y1), QPointF(x2, f(x2)));
        x1 = x2; 
	y1 = f(x2);
    }
}

// drawing Akima's or Splanes
void Window::DrawingApproximation(QPainter &painter, double a, double b, double dx,
                               int n, const double *X, const double *A, int m)
{
    double x1 = a;
	double y1 = EvaluationAkima(x1, a, b, n, X, A);
	double y2 = 0.0;
	if(m == 2) {
		y1 = EvaluationSplines(x1, a, b, n, X, A);
	}
    for (double x2 = x1 + dx; x2 <= b; x2 += dx) {
        y2 = EvaluationAkima(x2, a, b, n, X, A);
	if(m == 2) {
		y2 = EvaluationSplines(x2, a, b, n, X, A);
	}
        painter.drawLine(QPointF(x1, y1), QPointF(x2, y2));
        x1 = x2; 
	y1 = y2;
    }
}

// drawing error
void Window::DrawingError(QPainter &painter, double a, double b, double dx,
                       int n, const double *X, const double *A, int m)
{
    double x1 = a;
	double y1 = fabs(EvaluationAkima(x1, a, b, n, X, A) - f(x1));
	if(m == 2) {
		y1 = fabs(EvaluationSplines(x1, a, b, n, X, A) - f(x1));
	}
	double y2 = 0.0;
    for (double x2 = x1 + dx; x2 <= b; x2 += dx) {
        y2 = fabs(EvaluationAkima(x2, a, b, n, X, A) - f(x2));
	if(m == 2) {
		y2 = fabs(EvaluationSplines(x2, a, b, n, X, A)- f(x2));
	}
        painter.drawLine(QPointF(x1, y1), QPointF(x2, y2));
        x1 = x2; 
	y1 = y2;
    }
}


void Window::keyPressEvent(QKeyEvent *event)
{
    switch (event->key()) {
    case Qt::Key_0:
        k = (k + 1) % 7;
        func_id = k;
        change_func();
        break;       
    case Qt::Key_1:
        display = (display + 1) % 4;
        update();
        break;       
    case Qt::Key_2:
        scale++;
        update();
        break;      
    case Qt::Key_3:
        scale--;
        update();
        break;        
    case Qt::Key_4:
        n *= 2;
        if (n > 50){ n = 50; }
        update();
        break;      
    case Qt::Key_5:
        n /= 2;
        if (n < 2) { n = 2; }
        update();
        break;        
    case Qt::Key_6:
        perturbation++;
        update();
        break;        
    case Qt::Key_7:
        perturbation--;
        update();
        break;        
    default:
        QWidget::keyPressEvent(event);
    }
}

