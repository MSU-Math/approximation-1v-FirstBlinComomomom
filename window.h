#ifndef WINDOW_H
#define WINDOW_H

#include <QWidget>

class Window : public QWidget
{
    Q_OBJECT

  private:
    int func_id;
    const char *f_name;
    double a;
    double b;
    int n;
    int k;
    double (*f)(double);

    int display; // 7b
    int scale; // 7c
    int perturbation; // 7e
    double max_f;

    void DrawingFunction(QPainter &painter, double a, double b, double dx);
	void DrawingApproximation(QPainter &painter, double a, double b, double dx,
                               int n, const double *X, const double *A, int m);
	void DrawingError(QPainter &painter, double a, double b, double dx,
                       int n, const double *X, const double *A, int m);


  public:
    Window(QWidget *parent);

    QSize minimumSizeHint() const;
    QSize sizeHint() const;

    int parse_command_line(int argc, char *argv[]);

  public slots:
    void change_func();

  protected:
    void paintEvent(QPaintEvent *event);
    void keyPressEvent(QKeyEvent *event);
};

#endif
