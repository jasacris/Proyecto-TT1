#include "..\include\TimeUpdate.hpp"

Matrix TimeUpdate(Matrix &P, Matrix &Phi, Matrix &Qdt){
    return Phi*P*Phi.transpose() + Qdt;
}
