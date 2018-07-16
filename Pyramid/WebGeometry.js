// Перевод градусов в радианы
Math.radians = function(degrees) 
{
	return degrees * Math.PI / 180;
};

// Перевод радиан в градусы
Math.degrees = function(radians) 
{
	return radians * 180 / Math.PI;
};

// 2D - DIMENSION

// Функция-конструктор двумерных векторов Vector2D
function Vector2D(x, y)
{
	this[2];
	this[0] = x;
	this[1] = y;	
}

// Сложение двумерных векторов (this + vector)
// Функция возвращает новый вектор, получившийся в результат сложения
// с вектором this.
Vector2D.prototype.Add = function(vector) 
{
	var t1 = this[0] + vector[0];
	var t2 = this[1] + vector[1];
    return new Vector2D(t1, t2);
};

// Вычитание из вектора this вектора vector
// Функция возвращает новый вектор, получившийся в результат вычитания.
Vector2D.prototype.Subtract = function(vector) 
{
	var t1 = this[0] - vector[0];
	var t2 = this[1] - vector[1];
    return new Vector2D(t1, t2);
};

// Скалярное произведение векторов (this * vector)
// Функция возвращает результат произведения.
Vector2D.prototype.Dot = function(vector) 
{
	var out = this[0] * vector[0] + this[1] * vector[1];
    return out;
};

// Нормирование вектора this
Vector2D.prototype.Normer = function() 
{
	var norm = Math.sqrt(this[0] * this[0] + this[1] * this[1]);
    this[0] = this[0] / norm; 
	this[1] = this[1] / norm; 
};

// Вращение на плоскости вектора this на угол angle
// Функция возвращает новый вектор, получившийся в результат поворота.
Vector2D.prototype.Rotate = function(angle) 
{
	var mat = new Matrix2D();
	mat[0][0] = Math.cos(angle);
	mat[1][0] = -Math.sin(angle);
	// col 2
	mat[0][1] = Math.sin(angle);
	mat[1][1] = Math.cos(angle);	
	return mat.MultiplyMatrixVector(this);
};

// Перенос вектора this на вектор vector
// Фактически this = this + vector
Vector2D.prototype.Translate = function(vector)
{
	this[0] = this[0] + vector[0];
	this[1] = this[1] + vector[1];
}

// Matrix2D
// Функция-конструктор двумерных матриц Matrix2D
function Matrix2D(m00, m01, m10, m11)
{
	this[2];
	var row1 = [m00, m01];
	var row2 = [m10, m11];
	this[0] = row1;
	this[1] = row2;
}

// Функция возвращает величину определителя матрицы
Matrix2D.prototype.Det = function()
{
	if ( (this[0][0] == undefined) || (this[0][1] == undefined) ||
	   ( (this[1][0] == undefined) || (this[1][1] == undefined) ) )
	{
		return null;
	}
	var det = 0;
	det = this[0][0] * this[1][1] - this[0][1] * this[1][0];
	return det;
}

// Инверсия исходной матрицы this
Matrix2D.prototype.Inverse = function()
{
	var det;
	if ( (this[0][0] == undefined) || (this[0][1] == undefined) ||
	   ( (this[1][0] == undefined) || (this[1][1] == undefined) ) )
	{
		return null;
	}
	if((det = this.Det()) == 0.0) 
	{
		return null;
	}
	var outMat = new Matrix2D();
	outMat[0][0] =   this[0][0] / det;
	outMat[1][0] = - this[1][0] / det;
	outMat[0][1] = - this[0][1] / det;
	outMat[1][1] =   this[1][1] / det;
	return outMat;
}

// Функция возвращает транспонированную матрицу mat
Matrix2D.prototype.Transpose = function()
{
	var mat = new Matrix2D();
	var i, j;
	for(i = 0; i < 2; i++)
	{
		for( j = 0; j < 2; j++) 
		{
			mat[i][j] = this[j][i];
		}
	}
	return mat;	
}

// Функция возвращает вектор, получившийся в результате
// умножения матрицы this на вектор vector.
Matrix2D.prototype.MultiplyMatrixVector = function(vec)
{
	var t = [];
	t[0] = this[0][0]*vec[0] + this[0][1]*vec[1];
	t[1] = this[1][0]*vec[0] + this[1][1]*vec[1];
	var outVec = new Vector2D(t[0], t[1]);
	return outVec;
}

// Point2D

// Функция-конструктор точки на плоскости Point2D
// x и y задают координаты вновь создаваемой точки.
function Point2D(x, y)
{
	this[2];
	this[0] = x;
	this[1] = y;
}

// Функция возвращает точку с координатами равными сумме координт
// данной точки this и точки point
Point2D.prototype.Add = function(point) 
{
	var t1 = this[0] + point[0];
	var t2 = this[1] + point[1];
    return new Point2D(t1, t2);
};

//  Функция возвращает точку с координатами точки this 
// смещенными  на значение вектора vector
Point2D.prototype.Translate = function(vector)
{
	var v = new Vector2D(this[0], this[1]);
	v[0] = v[0] + vector[0];
	v[1] = v[1] + vector[1];
	var ptOut = new Point2D(v[0], v[1]);
	return ptOut;
}

// Функция возвращает значение расстояния между
// точками this и point
Point2D.prototype.Distance = function(point)
{
	var d0 = this[0] - point[0];
	var d1 = this[1] - point[1];
	return Math.sqrt(d0*d0 + d1*d1);
}

// Line2D
// Функция-конструктор прямой на плоскости Line2D
// Прямая задается точками point1 и point2
function Line2D(point1, point2)
{
	if (point1 == undefined)
	{
		this.directCos = [undefined, undefined];
		this.distOXY = undefined;
		return;
	}
	if (point2 == undefined)
	{
		this.directCos = [undefined, undefined];
		this.distOXY = undefined;
		return;
	}
	var dx, dy, norm, kon;
	dx = point2[0] - point1[0];
	dy = point2[1] - point1[1];

	kon = -point1[1]*dx + point1[0]*dy;

	norm = Math.sqrt(dx*dx + dy*dy);

	var cosx = dy/norm;
	var cosy = -dx/norm;
	var dist = kon/norm;
	this.directCos = [cosx, cosy];
	this.distOXY = dist;
}

// Функция создает прямую, перпендикулярную к исходной прямой.
Line2D.prototype.NormalLine = function()
{
	var line = new Line2D();
	line.directCos[0] = -this.directCos[1];
	line.directCos[1] =  this.directCos[0];
	line.distOXY = this.distOXY;
	return line;
}

// Функция создает прямую перпендикулярную к прямой this
// и проходящую через точку point.
// Данная функция НЕ является эквивалентом функции NormalLinе.
Line2D.prototype.CreateNormalLinePoint = function(point)
{
	var line = new Line2D();
	line.directCos[0] = -this.directCos[1];
	line.directCos[1] = this.directCos[0];
	line.distOXY = (point[0]) * (-this.directCos[1]) + (point[1]) * (this.directCos[0]);
	return line;
}
// Функция возвращает точку лежащую на пересечении прямых this и line.
Line2D.prototype.IntersectionTwoLines = function(line)
{
	var det, xp, yp;

	det = this.directCos[0]*line.directCos[1] - this.directCos[1]*line.directCos[0];

	if(det == 0.0) 
	{
		return null;
	}

	xp = this.distOXY * line.directCos[1] - line.distOXY * this.directCos[1];
	yp = this.directCos[0] * line.distOXY - line.directCos[0] * this.distOXY;

	var pt = new Point2D();
	pt[0] = xp / det;
	pt[1] = yp / det;

	return pt;
}

// Функция возвращает значение угла между прямыми this и line.
Line2D.prototype.Angle = function(line)
{
	var ang;
	var val = (this.directCos[0] * line.directCos[1] - this.directCos[1] * line.directCos[0]);
	var scal = (this.directCos[0] * line.directCos[0] + this.directCos[1] * line.directCos[1]);
	if(scal == 0.0)
	{
		ang = Math.PI/2;
	}
	else
	{
		ang = Math.atan2(val, scal);
	}
	return ang;
}

//  Данная прямая this создается по вектору направления vector 
// и точке point лежащей на прямой.
Line2D.prototype.CreateLineVectorPoint = function(vector, point)
{
	vector.Normer();
	this.directCos[0] = vector[1];
	this.directCos[1] = - vector[0];
	this.distOXY = vector[1]*point[0] - vector[0]*point[1];
}


// Функция возвращает значение расстояния прямой от заданной точки.
Line2D.prototype.Distance = function(point)
{
	return (this.directCos[0]*point[0] + this.directCos[1]*point[1] - this.distOXY);	
}

// Circle2D

// Функция-конструктор окружности на плоскости Circle2D
// Окружность задается своим центром point и величиной радиуса radius.
function Circle2D(point, radius)
{
	if ( (point == undefined) || (radius == undefined) )
	{
		this.center = [undefined, undefined];
		this.radius = undefined;
		return;
	}
	this.center = [point[0], point[1]];
	this.radius = radius;
}

// Функция определяет точки пересечения point1 и point2 двух окружностей.
// Одна окружность задается this, а вторая - circle2.
Circle2D.prototype.IntersectionTwoCircles = function(circle2, point1, point2)
{
	// point1, point2 - out cross points
	var midPt = new Point2D();
	var r1 = this.radius;
	var r2 = circle2.radius;
	
	var ptCenter1 = new Point2D(this.center[0], this.center[1]);
	var ptCenter2 = new Point2D(circle2.center[0], circle2.center[1]);
	
	// Distance between centers of circles
	var d = ptCenter1.Distance(ptCenter2);

	if ( d > (r1 + r2) ) 
	{
		return null;
	}

	if( ( d < r1) && ( (d + r2) < (r1 - 0.00001) ) ) // 0.00001 - очень маленькое значение ( ε )
	{
		return null;
	}

	if( (d < r2 ) && ( (d+r1) < (r2 - 0.00001) ) ) // 0.00001 - очень маленькое значение ( ε )
	{
		return null;
	}

	var a = (r1*r1 - r2*r2 + d*d) / (2*d);

	midPt[0] = ptCenter1[0] + a*(ptCenter2[0] - ptCenter1[0]) / d;
	midPt[1] = ptCenter1[1] + a*(ptCenter2[1] - ptCenter1[1]) / d;
	
	var t10 = ptCenter1[0];
	var t20 = ptCenter2[0];
	var tt = t10 + a*(t20 - t10)/d;

	if(Math.abs(r1) < Math.abs(a))
	{ 
		return null;
	}

	var h = Math.sqrt(r1*r1 - a*a);
	
	// Cross TwoLines 1
	point1[0] = midPt[0] - h*(ptCenter2[1] - ptCenter1[1])/d;
	point1[1] = midPt[1] + h*(ptCenter2[0] - ptCenter1[0])/d;

	// Cross TwoLines 2
	point2[0] = midPt[0] + h*(ptCenter2[1] - ptCenter1[1])/d;
	point2[1] = midPt[1] - h*(ptCenter2[0] - ptCenter1[0])/d;

	return 1;
}

//  Функция определяет точки пересечения point1 и point2 
// окружности this и прямой line.
Circle2D.prototype.IntersectionLineCircle = function(line, point1, point2)
{
	var pt1 = [];
	var pt2 = [];
	var tt = new Point2D(this.center[0], this.center[1]);
	if(line.Distance(tt) > this.radius)
	{
		return null;
	}

	var x0 = this.center[0];
	var y0 = this.center[1];
	var r  = this.radius;	
		
	var dist = line.distOXY; // !!!!!
	var norm = [line.directCos[0], line.directCos[1]];
	var x1 = dist * norm[0];
	var y1 = dist * norm[1];
	var pts = []; // Две точки пересечения
	
	if(Math.abs(norm[1]) > 0.00001)   // 0.00001 - очень маленькое значение ( ε )
	{
		var m = -norm[0]/norm[1];
		var a = 1 + m*m;
		var b = 2*(y1-y0-m*x1)*m-2*x0;
		var c = x0*x0+(y1-y0-m*x1)*(y1-y0-m*x1)-r*r;
		
		if(!QuadraticEquation(a, b, c, pts))
		{
			return null;
		}
		
		// Точка пересечения 1
		pt1[0] = pts[0];
		pt1[1] = y1 + m * (pts[0] - x1);
		point1[0] = pt1[0];
		point1[1] = pt1[1];
		
		// Точка пересечения 2
		pt2[0] = pts[1];
		pt2[1] = y1 + m * (pts[1] - x1);
		point2[0] = pt2[0];
		point2[1] = pt2[1];

		return 1;
	}
	else if(Math.abs(norm[0]) > 0.00001) // 0.00001 - очень маленькое значение ( ε )
	{
		var m = - norm[1]/norm[0];
		var a = 1 + m*m;
		var b = 2*(x1-x0-m*y1)*m-2*y0;
		var c = y0*y0+(x1-x0-m*y1)*(x1-x0-m*y1)-r*r;
		
		if(!QuadraticEquation(a, b, c, pts)) 
		{
			return null;
		}
		
		// The point of intersection 1
		pt1[0] = x1 + m * (pts[0] - y1);
		pt1[1] = pts[0];
		point1[0] = pt1[0];
		point1[1] = pt1[1];
		
		// The point of intersection 2
		pt2[0] = x1 + m * (pts[1] - y1);
		pt2[1] = pts[1];
		point2[0] = pt2[0];
		point2[1] = pt2[1];

		return 1;
	}
	return null;
}

// Решение квадратного уравнения.
// aa, bb, cc - коэффициенты и свободный член уравнения.
// rez - результат решения (массив из двух значений).
function QuadraticEquation(a, b, c, rez)
{
	var discr;
	discr = b*b - 4*a*c;
  
	if ( discr < 0.0 ) 
	{
		rez[0]= 0.0;
		rez[1]= 0.0;
		return null;
    }
	discr = Math.sqrt(discr);
    rez[0]= (-b + discr)/(2*a);
    rez[1]= (-b - discr)/(2*a);
	return 1;
}

// 3-DIMENSION

// Функция-конструктор векторов в пространстве Vector3D
function Vector3D(x, y, z)
{
	this[3];
	this[0] = x;
	this[1] = y;
	this[2] = z;	
}

// Сложение 3D векторов (this + vector)
// Функция возвращает новый вектор, получившийся в результат сложения.
Vector3D.prototype.Add = function(vector) 
{
	var t1 = this[0] + vector[0];
	var t2 = this[1] + vector[1];
	var t3 = this[2] + vector[2];
    return new Vector3D(t1, t2, t3);
};

// Вычитание из вектора this вектора vector
// Функция возвращает новый вектор, получившийся в результат вычитания.
Vector3D.prototype.Subtract = function(vector) 
{
	var t1 = this[0] - vector[0];
	var t2 = this[1] - vector[1];
	var t3 = this[2] - vector[2];
    return new Vector3D(t1, t2, t3);
};

// Скалярное произведение векторов (this * vector)
// Функция возвращает результат произведения.
Vector3D.prototype.Dot = function(vector) 
{
	var dot = this[0] * vector[0] + this[1] * vector[1] + this[2] * vector[2];
    return dot;
};

// Нормирование вектора this
Vector3D.prototype.Normer = function() 
{
	var norm = Math.sqrt(this[0] * this[0] + this[1] * this[1] + this[2] * this[2]);
    this[0] = this[0] / norm; 
	this[1] = this[1] / norm;
	this[2] = this[2] / norm;	
};

// Возвращает значение угла между векторами this и vector.
Vector3D.prototype.Angle = function(vector)
{
	var ang;
	var scal = this[0] * vector[0] + this[1] * vector[1] + this[2] * vector[2];
	var norm1 = Math.sqrt(this[0] * this[0] + this[1] * this[1] + this[2] * this[2]);
	var norm2 = Math.sqrt(vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2]);
	
	scal = scal / (norm1 * norm2); 
       
	if(scal > 1.0)
	{
		scal = 1.0;
	}
	if (scal <= -1.0)
	{
		scal = -1.0;
	}
	
	ang = Math.acos(scal);
	return ang;
}

// Возвращает вектор, являющийся векторным произведением
// векторов this и vector.
Vector3D.prototype.Cross = function(vector)
{
	var vec = new Vector3D();
	vec[0] = this[1] * vector[2] - this[2] * vector[1];
	vec[1] = this[2] * vector[0] - this[0] * vector[2];
	vec[2] = this[0] * vector[1] - this[1] * vector[0];
	return vec;
}

// Вектор this поворачивается на величину определяемую матрицей rot.
Vector3D.prototype.Rotate = function(rot)
{
	var coord = new Vector3D(this[0], this[1], this[2]);
	var vec = rot.MultiplyMatrixVector(coord);	
	//this[0] = result[0];
	//this[1] = result[1];
	//this[2] = result[2];
	return vec;
}

// Matrix3D
// Функция-конструктор трехмерных матриц Matrix3D
function Matrix3D(m0x, m0y, m0z, m1x, m1y, m1z, m2x, m2y, m2z)
{
	// m0x  m1x  m2x
	// m0y  m1y  m2y
	// m0z  m1z  m2z
	var col0 = [m0x, m0y, m0z];
	var col1 = [m1x, m1y, m1z];
	var col2 = [m2x, m2y, m2z];
	this[0] = col0;
	this[1] = col1;
	this[2] = col2;
}

// Функция возвращает величину определителя матрицы
Matrix3D.prototype.Det = function()
{
	if ( (this[0][0] == undefined) || (this[0][1] == undefined) || (this[0][2] == undefined) ||
		 (this[1][0] == undefined) || (this[1][1] == undefined) || (this[1][2] == undefined) ||
		 (this[2][0] == undefined) || (this[2][1] == undefined) || (this[2][2] == undefined) )
	{
		return null;
	}
	var det = 0;
	det  = this[0][0] * this[1][1] * this[2][2];
	det += this[0][1] * this[1][2] * this[2][0];
	det += this[1][0] * this[2][1] * this[0][2];
	det -= this[2][0] * this[1][1] * this[0][2];
	det -= this[0][1] * this[1][0] * this[2][2];
	det -= this[0][0] * this[2][1] * this[1][2];

	return det;
}

// Инверсия исходной матрицы this
// Функция возвращает обратную матрицу к матрице this
Matrix3D.prototype.Inverse = function()
{
	var det, tmp;
	var i, j, m, n, r, s;
	det = this.Det();
	if((det == 0) || (det == null))
	{
		return null;
	}	
	var mat = new Matrix3D();
	for(i = 0; i < 3; i++)
	{
		m = (i+2) % 3;   // m = i-1;
		n = (i+1) % 3;	// n = i+1;

		for(j = 0; j < 3; j++)
		{
			r = (j+2) % 3;	// r = j-1;
			s = (j+1) % 3;	// s = j+1;
				
			tmp = this[m][r]*this[n][s] - this[m][s]*this[n][r];
			mat[j][i] = tmp / det;
		}
	}
	return mat;
}

// Возвращает вектор, полученный путем умножения матрицы this на вектор vector.
// vecOut = this * vector;
Matrix3D.prototype.MultiplyMatrixVector = function(vec)
{
	var vecOut = new Vector3D();
	vecOut[0] = this[0][0]*vec[0] + this[0][1]*vec[1] + this[0][2]*vec[2];
	vecOut[1] = this[1][0]*vec[0] + this[1][1]*vec[1] + this[1][2]*vec[2];
	vecOut[2] = this[2][0]*vec[0] + this[2][1]*vec[1] + this[2][2]*vec[2];
	return vecOut;
}

// Возвращает матрицу, полученную путем умножения матрицы this на матрицу mat.
// matOut = this * mat
Matrix3D.prototype.MultiplyTwoMatrices = function(mat) 
{
	var matOut = new Matrix3D();
	var i, j, k;
	var sum;

	for(i = 0; i < 3; i++)
	{
		for(j = 0; j < 3; j++)
		{
			sum = 0.0;
			for(k = 0; k < 3; k++)
			{
				sum += this[i][k]*mat[k][j];
			}
			matOut[i][j] = sum;
		}
	}
	return matOut;
}

// Возвращает транспонированную матрицу mat.
Matrix3D.prototype.Transpose = function(mat)
{
	var matOut = new Matrix3D();
	var i, j;
	for(i = 0; i < 3; i++)
	{
		for( j = 0; j < 3; j++) 
		{
			matOut[i][j] = mat[j][i];
		}
	}
	return matOut;
}


// Устанавливает значения матрицы по столбцам
// col1, col2, col3 - 3D-векторы
Matrix3D.prototype.MatrixSetVectorInColumn = function(col1, col2, col3)
{
	var i;
	for(i = 0; i < 3; i++)
	{
	  this[i][0] = col1[i]; 
	  this[i][1] = col2[i];
	  this[i][2] = col3[i];
	}
	return this;
}

// Устанавливает значения матрицы по рядам
// vec0, vec1, vec2 - 3D-векторы
Matrix3D.prototype.MatrixSetVectorInRow = function(vec0, vec1, vec2)
{
	var col0 = [vec0[0], vec0[1], vec0[2]];
	var col1 = [vec1[0], vec1[1], vec1[2]];
	var col2 = [vec2[0], vec2[1], vec2[2]];
	
	this[0] = col0;
	this[1] = col1;
	this[2] = col2;
	
	return this;
}

// Единичная матрица
Matrix3D.prototype.UnitaryMatrix = function()
{ 
	var oneX = new Vector3D(1, 0, 0);
	var oneY = new Vector3D(0, 1, 0);
	var oneZ = new Vector3D(0, 0, 1);
	this.MatrixSetVectorInColumn(oneX, oneY, oneZ);
	return this;
}

// Функция возвращает инвертированную матрицу this умноженную
// на вектор vector.
// Данная функция используется в CrossThreePlanes.
Matrix3D.prototype.InvertMatrixMultiply = function (vector)
{
	var vecOut = new Vector3D();
	if(this.Det() == 0.0)
	{
		return null;
	}
	
	var matrix = this.Inverse();
	vecOut = matrix.MultiplyMatrixVector(vector);
	return vecOut;
}

// Матрица поворота вокруг оси X на угол angle
Matrix3D.prototype.RotX = function(angle)
{
	// col 1
	this[0][0] = 1.0;
	this[1][0] = 0.0;
	this[2][0] = 0.0;
	// col 2
	this[0][1] = 0.0;
	this[1][1] = Math.cos(angle);
	this[2][1] = Math.sin(angle);
	// col 3
	this[0][2] = 0.0;
	this[1][2] = - Math.sin(angle);
	this[2][2] = Math.cos(angle);
	return this;
}

// Матрица поворота вокруг оси Y на угол angle
Matrix3D.prototype.RotY = function(angle)
{
	//col 1
	this[0][0] = Math.cos(angle);
	this[1][0] = 0.0;
	this[2][0] = - Math.sin(angle);
	// col 2
	this[0][1] = 0.0;
	this[1][1] = 1.0;
	this[2][1] = 0.0;
	// col 3
	this[0][2] = Math.sin(angle);
	this[1][2] = 0.0;
	this[2][2] = Math.cos(angle);
	return this;
}

// Матрица поворота вокруг оси Z на угол angle
Matrix3D.prototype.RotZ = function(angle)
{
	// col 1
	this[0][0] = Math.cos(angle);
	this[1][0] = Math.sin(angle);
	this[2][0] = 0.0;
	// col 2
	this[0][1] = - Math.sin(angle);
	this[1][1] = Math.cos(angle);
	this[2][1] = 0.0;
	// col 3
	this[0][2] = 0.0;
	this[1][2] = 0.0;
	this[2][2] = 1.0;
	return this;
}

// Установка компонент матрицы вращения (установка компонент матрицы this)
// по трем векторам vec1, vec2 и vec3.
Matrix3D.prototype.SetRotateMatrix = function(vec1, vec2, vec3)
{
	vec1.Normer();
	vec2.Normer();
	vec3.Normer();
	
	//var rot = new Matrix3D();
	this.MatrixSetVectorInColumn(vec1, vec2, vec3);
	var mat = new Matrix3D();
	mat = this.Transpose(this);
	
	var i, j;
	for(i = 0; i < 3; i++)
	{
		for( j = 0; j < 3; j++) 
		{
			this[i][j] = mat[i][j];
		}
	}
	return mat;
}

// Point3D
// Функция-конструктор 3D-точки Point3D
// x, y, z - координаты точки.
function Point3D(x, y, z)
{
	this[3];
	this[0] = x;
	this[1] = y;
	this[2] = z;
}

// Возвращает точку, являющюся суммой двух точек this и point.
Point3D.prototype.Add = function(point) 
{
	var t1 = this[0] + point[0];
	var t2 = this[1] + point[1];
	var t3 = this[2] + point[2];
    return new Point3D(t1, t2, t3);
};

// Возвращает точку this перенесенную на вектор vector
Point3D.prototype.Translate = function(vector)
{
	var v = new Vector3D(this[0], this[1], this[2]);
	v[0] = v[0] + vector[0];
	v[1] = v[1] + vector[1];
	v[2] = v[2] + vector[2];
	var ptOut = new Point3D(v[0], v[1], v[2]);
	return ptOut;
}

// Возвращает расстояние между точкой this и точкой point.
Point3D.prototype.Distance = function(point)
{
	var d0 = this[0] - point[0];
	var d1 = this[1] - point[1];
	var d2 = this[2] - point[2];
	return Math.sqrt(d0*d0 + d1*d1 + d2*d2);
}

// Функция возвращает точку this, повернутую матрицей rot.
Point3D.prototype.Rotate = function(rot)
{
	var vec = new Vector3D(this[0], this[1], this[2]);
	var vecOut = rot.MultiplyMatrixVector(vec)
	var ptOut = new Point3D(vecOut[0], vecOut[1], vecOut[2]);
	return ptOut;
}

// Функция возвращает площадь треугольника заданного точками pt1, pt2, и pt3
function SquareTriangle(pt1, pt2, pt3)
{
	var mat1 = new Matrix3D();
	var mat2 = new Matrix3D();
	var mat3 = new Matrix3D();
	
	var det1, det2, det3;
	
	// 1 matrix
	mat1[0][0] = pt1[1];
	mat1[1][0] = pt2[1];
	mat1[2][0] = pt3[1];
	mat1[0][1] = pt1[2];
	mat1[1][1] = pt2[2];
	mat1[2][1] = pt3[2];
	mat1[0][2] = 1.0;
	mat1[1][2] = 1.0;
	mat1[2][2] = 1.0;
	det1 = mat1.Det();
	// 2 matrix
	mat2[0][0] = pt1[2];
	mat2[1][0] = pt2[2];
	mat2[2][0] = pt3[2];
	mat2[0][1] = pt1[0];
	mat2[1][1] = pt2[0];
	mat2[2][1] = pt3[0];
	mat2[0][2] = 1.0;
	mat2[1][2] = 1.0;
	mat2[2][2] = 1.0;
	det2 = mat2.Det();
	// 3 matrix
	mat3[0][0] = pt1[0];
	mat3[1][0] = pt2[0];
	mat3[2][0] = pt3[0];
	mat3[0][1] = pt1[1];
	mat3[1][1] = pt2[1];
	mat3[2][1] = pt3[1];
	mat3[0][2] = 1.0;
	mat3[1][2] = 1.0;
	mat3[2][2] = 1.0;
	det3 = mat3.Det();

	return Math.sqrt(det1*det1 + det2*det2 + det3*det3) / 2.0;
}

// Line3D
// функция создает новую прямую в пространстве по двум точкам pt1 и pt2.
function Line3D(pt1, pt2)
{
	this.directCos = [3];
	if ( (pt1 == undefined) || (pt2 == undefined) )
	{

		this.directCos[0] = undefined;
		this.directCos[1] = undefined;
		this.directCos[2] = undefined;
		var pt = new Point3D();
		this.pointOnLine = pt;
		return;
	}
	var t = new Vector3D();
	t[0] = pt2[0] - pt1[0];
	t[1] = pt2[1] - pt1[1];
	t[2] = pt2[2] - pt1[2];
	t.Normer();
	this.directCos[0] = t[0];
	this.directCos[1] = t[1];
	this.directCos[2] = t[2];

	this.pointOnLine = pt1;	
}	

// функция создает прямую в пространстве направленную вдоль 
// вектора vec и проходящую через точку pt.
Line3D.prototype.CreateLineVectorPoint = function(vec, pt)
{
	vec.Normer();
	this.directCos = [vec[0], vec[1], vec[2]];
	this.pointOnLine = pt;
	return this;
}

// Функция вращает прямую this в пространстве посредством
// матрицы rot
Line3D.prototype.Rotate = function(rot)
{
	var t = new Vector3D();
	t[0] = this.directCos[0];
	t[1] = this.directCos[1];
	t[2] = this.directCos[2];
	var vecRot = t.Rotate(rot);
	this.directCos[0] = vecRot[0];
	this.directCos[1] = vecRot[1];
	this.directCos[2] = vecRot[2];
	var pt = this.pointOnLine.Rotate(rot);
	this.pointOnLine = pt;
}

// Функция возвращает угол между векторами определяющими
// направления прямых this и line.
Line3D.prototype.Angle = function(line)
{
	var ang;
	ang = Math.acos(this.directCos[0]*line.directCos[0] 
	           + this.directCos[1]*line.directCos[1] 
			   + this.directCos[2]*line.directCos[2]);
	return ang;
}

// Функция-конструктор трехмерной плоскости.
// Плоскость задается вектором нормали к плоскости norm
// и расстоянием dist от начала координат.
function Plane3D(norm, dist)
{
	this.directCos = [3];
	if ( (norm == undefined) || (dist == undefined) )
	{
		this.distOXYZ = undefined;
		this.directCos[0] = undefined;
		this.directCos[1] = undefined;
		this.directCos[2] = undefined;	
		return;
	}
	this.directCos[0] = norm[0];
	this.directCos[1] = norm[1];
	this.directCos[2] = norm[2];
	this.distOXYZ = dist;
}

// Функция создает плоскость по направляющему вектору norm и расстоянию от начала
// координат. Делает тоже, что и конструктор плоскости.
Plane3D.prototype.CreatePlaneNormalDistOXYZ = function(norm, dist)  
{
	norm.Normer();
	this.directCos[0] = norm[0];
	this.directCos[1] = norm[1];
	this.directCos[2] = norm[2];
	this.distOXYZ = dist;
}

// Нормировка уравнения плоскости.
Plane3D.prototype.Normer = function() 
{
	var norm;

	norm = Math.sqrt(this.directCos[0]*this.directCos[0] + this.directCos[1]*this.directCos[1] + 
								this.directCos[2]*this.directCos[2] );
	
	if(norm == 0.0) 
		return null;
	
	this.directCos[0] = this.directCos[0] / norm;
	this.directCos[1] = this.directCos[1] / norm;
	this.directCos[2] = this.directCos[2] / norm;
	this.distOXYZ = this.distOXYZ / norm;

	return 1;
}

// Функция создает плоскость проходящую через три точки pt1, pt2 и pt3.
Plane3D.prototype.CreatePlaneThreePoints = function(pt1, pt2, pt3)
{
	var col1 = new Vector3D(pt1[0], pt2[0], pt3[0]);
	var col2 = new Vector3D(pt1[1], pt2[1], pt3[1]);
	var col3 = new Vector3D(pt1[2], pt2[2], pt3[2]);
	var col4 = new Vector3D(1,1,1);
	var mat = new Matrix3D();

	var v = new Vector3D() ;
	
	// Уравнение Ax + By + Cz = D

	// Расстояние D от (0, 0, 0)
	mat.MatrixSetVectorInColumn(col1, col2, col3);
	this.distOXYZ = mat.Det();

	// A
	mat.MatrixSetVectorInColumn(col2, col3, col4);
	v[0] = mat.Det();
	// B
	mat.MatrixSetVectorInColumn(col1, col3, col4);
	v[1] = -mat.Det();
	// C
	mat.MatrixSetVectorInColumn(col1, col2, col4);
	v[2] = mat.Det();

	this.directCos[0] = v[0];
	this.directCos[1] = v[1];
	this.directCos[2] = v[2];
	
	this.Normer();
}

// Создание плоскости по ее нормальному вектору и точке, принадлежащей плоскости.
Plane3D.prototype.CreatePlaneNormalVectorPoint = function(vec, pt)
{
	vec.Normer();
	this.directCos = vec;
	this.distOXYZ = pt[0]*vec[0] + pt[1]*vec[1] + pt[2]*vec[2];
	this.Normer();
}

// Создание плоскости по вектору vec1, параллельному плоскости, 
// и двум точкам pt1 и pt2, принадлежащим плоскости.
Plane3D.prototype.CreatePlaneVectorTwoPoints = function(vec1, pt1, pt2)
{
	var vec2 = new Vector3D();
	var normal = new Vector3D();
	vec1.Normer();
	vec2 = (new Vector3D(pt2[0], pt2[1], pt2[2])).Subtract(new Vector3D(pt1[0], pt1[1], pt1[2]));
	normal = vec1.Cross(vec2);
	normal.Normer();
	this.CreatePlaneNormalVectorPoint(normal, pt1);
	this.Normer();
	return this;
}

// Возвращает точку пересечения плоскости this с плоскостями pl2 и pl3.
Plane3D.prototype.IntersectionThreePlanes = function(pl2, pl3)
{
	var pl1 = new Plane3D();
	pl1.directCos[0] = this.directCos[0];
	pl1.directCos[1] = this.directCos[1];
	pl1.directCos[2] = this.directCos[2];
	pl1.distOXYZ = this.distOXYZ;

	var mat = new Matrix3D();
	var i;
	for(i = 0; i < 3; i++)
	{
		mat[0][i] = pl1.directCos[i];
		mat[1][i] = pl2.directCos[i];
		mat[2][i] = pl3.directCos[i];
	}

	var vec = new Vector3D( pl1.distOXYZ, pl2.distOXYZ, pl3.distOXYZ );
	var vecOut = mat.InvertMatrixMultiply(vec);
	var pt = new Point3D(vecOut[0], vecOut[1], vecOut[2]);
	return pt;
}

// Функция возвращает вектор, имеющий направление прямой по которой
// пересекаются плоскости this и pl2 
Plane3D.prototype.VectorIntersectionTwoPlanes = function(pl2)
{
	var pl1 = new Plane3D();
	pl1.directCos[0] = this.directCos[0];
	pl1.directCos[1] = this.directCos[1];
	pl1.directCos[2] = this.directCos[2];
	pl1.distOXYZ = this.distOXYZ

	var vec = new Vector3D();
	vec[0] =   pl1.directCos[1]*pl2.directCos[2] 
	         - pl1.directCos[2]*pl2.directCos[1];
	vec[1] =   pl1.directCos[2]*pl2.directCos[0] 
	         - pl1.directCos[0]*pl2.directCos[2];
	vec[2] =   pl1.directCos[0]*pl2.directCos[1] 
	         - pl1.directCos[1]*pl2.directCos[0];
	vec.Normer();
	return vec;	
}

// Функция при помощи матрицы rot вращает плоскость this
Plane3D.prototype.Rotate = function(rot)
{
	var vecRotate = new Vector3D(this.directCos[0], this.directCos[1], this.directCos[2]);
	this.directCos = rot.MultiplyMatrixVector(vecRotate);	
}

// Функция возвращает расстояние от точки pt до плоскости this.
// return distance point from plane
Plane3D.prototype.DistancePoint = function(pt)
{
	return this.directCos[0]*pt[0] + this.directCos[1]*pt[1] + this.directCos[2]*pt[2] - this.distOXYZ;	
}

// Функция возвращает угол между плоскостью this и плоскостью OXY. 
// return angle between current plane and plane OXY
Plane3D.prototype.Slope = function()
{
	var d = Math.sqrt(this.directCos[0]*this.directCos[0] + this.directCos[1]*this.directCos[1]);
	var ang = Math.atan2(d, this.directCos[2]);
	if(ang > Math.PI/2) 
	{
		ang = ang - Math.PI
	}
	return ang;
}

// Функция возвращает угол между плоскостью this и плоскостью pl.
// return angle between two planes
Plane3D.prototype.Angle = function(pl)
{
	var ang;
	ang = Math.acos(this.directCos[0]*pl.directCos[0] 
	             + this.directCos[1]*pl.directCos[1]
	      	     + this.directCos[2]*pl.directCos[2]);
	return ang;
}

// Функция возвращает вектор перпендикулярный к плоскости this.
Plane3D.prototype.Normal = function()
{
	return new Vector3D(this.directCos[0], this.directCos[1], this.directCos[2]);
}

// Функция возвращает расстояние от плоскости this до начала координат.
Plane3D.prototype.Distance = function()
{
	return this.distOXYZ;
}

// Функция возвращает точку пересечения прямой this с плоскостью pl.
Line3D.prototype.IntersectionLinePlane = function(pl)
{
	//  x = x1 + a0 * t; 

	//var norm_plane = new Vector3D(pl.normal[0], pl.normal[1], pl.normal[2]);
	var normal = pl.Normal();
	var pt2 = new Vector3D(this.pointOnLine[0], this.pointOnLine[1], this.pointOnLine[2]);

	var val = pl.Distance() - (pt2.Dot(normal));
	var vec = new Vector3D(this.directCos[0], this.directCos[1], this.directCos[2] );
	//var scal = this.m_directCos.dot(norm_plane);
	var scal = vec.Dot(normal);

	if(scal == 0) 
	{
		return null; // Прямая лежит на плоскости или параллельна ей
	}

	var t = val / scal;

	var pt1 = new Point3D();
	pt1[0] = pt2[0] + this.directCos[0] * t; 
	pt1[1] = pt2[1] + this.directCos[1] * t; 
	pt1[2] = pt2[2] + this.directCos[2] * t; 

	return pt1;
}

// Предполагается, что точки pt1 и pt2 имеют одинаковую величину
// координаты z (иными словами pt1[2] = pt2[2]).
// Рассчитывается двумерный вектор norm2d и затем определяется 
// возвращаемая плоскость имеющая наклон angle к плоскости OXY.
function Facet(angle, pt1, pt2, pt3)	
{
	var norm2d = new Vector2D(pt1[1] - pt2[1], pt2[0] - pt1[0]);
	norm2d.Normer();
	var x = Math.sin(angle) * norm2d[0];
	var y = Math.sin(angle) * norm2d[1];
	var z = Math.cos(angle);
	var normPlaneVector = new Vector3D(x, y, z);
	var plane = new Plane3D();
	plane.CreatePlaneNormalVectorPoint(normPlaneVector, pt3);
	return plane;
} 

//  Данная функция создает плоскость наклоненную относительно плоскости OXY
// на угол angle.
// Точки pt1 и pt2 должны иметь одинаковое значение координаты z
//   (иными словами pt1[2] = pt2[2]).
// Рассчитывается двумерный вектор norm2d и затем определяется 
// плоскость имеющая наклон angle к плоскости OXY.
// Данная функция аналогична функции Facet.
Plane3D.prototype.CreateInclinePlane = function(angle, pt1, pt2, pt3)
{
	var norm2d = new Vector2D(pt1[1] - pt2[1], pt2[0] - pt1[0]);
	norm2d.Normer();
	var x = Math.sin(angle) * norm2d[0];
	var y = Math.sin(angle) * norm2d[1];
	var z = Math.cos(angle);
	var normPlaneVector = new Vector3D(x, y, z);
	normPlaneVector.Normer();
	this.CreatePlaneNormalVectorPoint(normPlaneVector, pt3);
	this.Normer();
}
