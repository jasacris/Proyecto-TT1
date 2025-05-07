#include "..\include\TimeUpdate.hpp"

Matrix& TimeUpdate(Matrix &P, Matrix &Phi, double Qdt){
    return Phi*P*transponse(Phi) + Qdt;
}