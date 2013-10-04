<?php

namespace Sole\Model;

class Sole {

    /**
     * Система линейных алгебраических уравнений
     */
    private $_system;

    /**
     * @param array $system СЛАУ
     * @return Sole
     * Установить систему уравнений
     */
    public function setSystem($system) {
        $this->_system = $system;
        return $this;
    }

    /**
     * @return array СЛАУ
     * Получить систему уравнений
     */
    public function getSystem() {
        return $this->_system;
    }

    /**
     * Получить размерность системы уравнений
     * @return integer размерность системы уравнений
     */
    public function getDimension() {
        return count($this->_system);
    }

    /**
     * Проверить систему уравнений на корректность заполнения
     * @return boolean корректность заполнения системы
     */
    public function isSystemValid() {
        $dimension = $this->getDimension();
        if ($dimension < 1)
            return FALSE;
        $expectedCount = $dimension + 1; //размерность +1 колонка результата
        foreach ($this->_system as $equation)
            if (count($equation) != $expectedCount)
                return FALSE;
        return TRUE;
    }

    /**
     * Получить  СЛАУ, приведенную к треугольному виду
     * @return array СЛАУ, приведенная к треугольному виду
     */
    public function getTriangularSystem() {
        if (!$this->isSystemValid())
            throw new \Exception("СЛАУ не прошла валидацию");
        $system = $this->_system;
        $dimension = $this->getDimension();
        for ($i = 0; $i < $dimension; $i++) {
            //var_dump($system);
            $coefficient = $system[$i][$i]; //коэфицент для деления из текущего уравнения
            $newEquation = array();
            foreach ($system[$i] as $item)
                $newEquation[] = $item / $coefficient; //новое уравнение все элементы /коэфицент
            $system[$i] = $newEquation; // замена i-го уравнения в системе
            for ($j = $i + 1; $j < $dimension; $j++) {
                $foreignCofficient = $system[$j][$i]; //коэфицент для исключения из уравнения
                for ($k = $i; $k < $dimension + 1; $k++) {
                    $system[$j][$k] = $newEquation[$k] * $foreignCofficient + $system[$j][$k] * -1;
                }
            }
        }
        return $system;
    }

    /**
     * Получить  решение СЛАУ методом Гаусса
     * @return array решение СЛАУ методом Гаусса
     */
    public function getSolutionGaussMethod() {
        $triangularSystem = $this->getTriangularSystem();
        $dimension = $this->getDimension();
        $solution = array();
        for ($i = $dimension - 1; $i >= 0; $i--) {
            $sol = $triangularSystem[$i][$dimension];
            for ($j = $i; $j < $dimension - 1; $j++) {
                $sol = $sol - $triangularSystem[$i][$j + 1] * $solution[$j + 1];
            }
            $solution[$i] = $sol;
        }
        return array_reverse($solution);
    }

    /**
     * Получить  СЛАУ, приведенную к симметричному виду
     * @return array СЛАУ, приведенная к симметричному виду
     */
    public function getSymmetricSystem() {
        $dimension = $this->getDimension();
        $c = array();
        for ($i = 1; $i <= $dimension; $i++) {
            for ($j = 1; $j <= $dimension; $j++) {
                $system[$i][$j] = $this->_system[$i - 1][$j - 1]; //Переход от си стандарта нумерации к натуральной
            }
        }//Умножение транспонированной на матрицу
        for ($i = 1; $i <= $dimension; $i++) {
            for ($j = 1; $j <= $dimension; $j++) {
                $summ = 0;
                for ($t = 1; $t <= $dimension; $t++) {
                    $summ = $system[$t][$j] * $system[$t][$i] + $summ;
                }
                $c[$j][$i] = $summ;
            }
        }//Умножение транспонированной на вектор свободных членов
        for ($i = 1; $i <= $dimension; $i++) {
            $y[$i] = 0;
            for ($j = 1; $j <= $dimension; $j++) {
                $y[$i]+=$system[$j][$i] * $this->_system[$j - 1][$dimension];
            }
        }
        $c['y'] = $y;
        return $c;
    }

    /**
     * Получить  решение СЛАУ методом Квадратных корней
     * @return array решение СЛАУ методом Квадратных корней
     */
    public function getSolutionCholeskyMethod() {
        if (!$this->isSystemValid())
            throw new \Exception("СЛАУ не прошла валидацию");
        $dimension = $this->getDimension();
        for ($i = 1; $i <= $dimension; $i++) {
            for ($j = 1; $j <= $dimension; $j++) {
                $l[$i][$j] = 0; //заполнение незначащих элементов матриц разложений нулями(для наглядности view)
                $lt[$i][$j] = 0;
                $system[$i][$j] = $this->_system[$i - 1][$j - 1];
            }
        }
        $system = $this->getSymmetricSystem();
        for ($i = 1; $i <= $dimension; $i++) {
            if ($i == 1)
                $uii = sqrt($system[$i][$i]);
            else {
                $sum = 0;
                for ($k = 1; $k < $dimension; $k++) {
                    $sum+=$l[$k][$i] * $l[$k][$i];
                }
                $uii = sqrt($system[$i][$i] - $sum); //Вычисление U(i,i) элемента через сумму
                $sqrt = $system[$i][$i] - $sum;
            }
            $l[$i][$i] = $uii;
            $lt[$i][$i] = $uii;
            for ($j = 2; $j <= $dimension; $j++) {
                $sum = 0;
                for ($k = 1; $k <= $i - 1; $k++) {
                    $sum+=$l[$k][$i] * $l[$k][$j];
                }
                $uij = ($system[$i][$j] - $sum) / $uii; //Вычисление U(i,j) элемента через сумму и U(i,i)
                $l[$i][$j] = $uij; //запись в матрицы разложений
                $lt[$j][$i] = $uij;
            }
        }
        $b = $system['y'];
        for ($i = 1; $i <= $dimension; $i++) {
            $mid = $b[$i];
            for ($j = $i; $j > 1; $j--) {
                $mid = $mid - $y[$j - 1] * $l[$j - 1][$i];
            }
            $y[$i] = $mid / $l[$i][$i]; // y - вектор решения U*y=b
        }
        for ($i = $dimension; $i > 0; $i--) {
            $mid = $y[$i];
            for ($j = $i; $j < $dimension; $j++) {
                $mid = $mid - $x[$j + 1] * $lt[$j + 1][$i];
            }
            $x[$i] = $mid / $lt[$i][$i]; // x - вектор решения Ut*x=y
        }
        $x = array_reverse($x);
        $x['l'] = $l; //Для полноты view
        $x['lt'] = $lt;
        $x['y'] = $y;
        return $x;
    }

}
