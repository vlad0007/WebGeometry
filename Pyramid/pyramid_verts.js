// Расчет координат вершин огранки (модели).
function VerticesCalculation()
{
	InitGirdle ();
	
	// Вспомогательные массивы
	var crown = [4];
	var pavil = [1];
	
	// Единичный вектор расположенный по оси OZ.
	var Z1 = new Vector3D(0, 0, 1);
	
	//    Создание плоскостей короны A, B, C и D.
	// Важна последовательность задания координат точек для
    // параметров pt1 и pt2 в функциях  CreateInclinePlane и Facet.
	// Принципиальное значение имеет то, какая точка идет первой, 
	// а какая - второй. От этого зависит вперед или назад (влево или вправо)
	// повернется на угол angleA (или angleB) соответствующая плоскость.
	// Плоскости можно создавать как функцией CreateInclinePlane 
	// так и функцией Facet.
	
	// Создаем плоскрость в которой расположена грань короны A.
	var A = new Plane3D();	
	A.CreateInclinePlane(angleA, girdle[3], girdle[0], girdle[0]);
	// Создаем плоскрость в которой расположена грань короны B.									   
	var B = Facet(angleB, girdle[0], girdle[1], girdle[1]);
	// Создаем плоскрость в которой расположена грань короны C.
	var C = new Plane3D();	
	C.CreateInclinePlane(angleA, girdle[1], girdle[2],  girdle[2]);
	// Создаем плоскрость в которой расположена грань короны D.									   
	var D = Facet(angleB, girdle[2], girdle[3], girdle[3]);

	// Создаем плоскость в которой лежит площадка Table.
	var Table = new Plane3D();
	Table.CreatePlaneNormalDistOXYZ(Z1, hCrown + r/2);
	
	// Координаты вершин короны 0, 1, 2 и 3 находим как точки
	// пересечения соответствующих граней короны с плоскостью
	// в которой лежит площадка Table.
	crown[0] = Table.IntersectionThreePlanes(A, B);
	crown[1] = Table.IntersectionThreePlanes(B, C);
	crown[2] = Table.IntersectionThreePlanes(C, D);
	crown[3] = Table.IntersectionThreePlanes(D, A);
/*	
	//  После того как были найдены координаты вершины crown[0]
	// координаты вершин crown[1], crown[2] и crown[3]
	// можно было найти из соображений симметрии огранки.
	crown[1] = new Point3D(  crown[0][0], - crown[0][1],   crown[0][2]);
	crown[2] = new Point3D(- crown[0][0], - crown[0][1],   crown[0][2]);
	crown[3] = new Point3D(- crown[0][0],   crown[0][1],   crown[0][2]);
*/
	// Создаем плоскрость в которой расположена грань павильона PavA.
	var PavA = Facet(-anglePav, girdle[7], girdle[4], girdle[4]);
								
	// Создаем прямую проходящую вертикально через начало координат.
	var line = new Line3D(new Point3D(0, 0, 0), new Point3D(0, 0, 1));
	
	// Единственную в данной огранке вершину павильона pavil[0]
	// находим как точку пересечения прямой line с плоскостью 
	// в которой лежит грань PavA.
	pavil[0] = line.IntersectionLinePlane(PavA);	
	
	// Заполняем массив vertices.
	var i;
	for(i = 0; i < 4; i++)
	{
		vertices.push(crown[i][0]);
		vertices.push(crown[i][1]);
		vertices.push(crown[i][2]);
	}
	for(i = 0; i < 8; i++)
	{
		vertices.push(girdle[i][0]);
		vertices.push(girdle[i][1]);
		vertices.push(girdle[i][2]);
	}
	for(i = 0; i < 1; i++)
	{
		vertices.push(pavil[i][0]);
		vertices.push(pavil[i][1]);
		vertices.push(pavil[i][2]);
	}								 
}

function InitGirdle ()
{
	girdle[0] = new Point3D(  lw * 0.5,   0.5,   r/2);
	girdle[1] = new Point3D(  lw * 0.5, - 0.5,   r/2);
	girdle[2] = new Point3D(- lw * 0.5, - 0.5,   r/2);
	girdle[3] = new Point3D(- lw * 0.5,   0.5,   r/2);
	girdle[4] = new Point3D(  lw * 0.5,   0.5, - r/2);
	girdle[5] = new Point3D(  lw * 0.5, - 0.5, - r/2);
	girdle[6] = new Point3D(- lw * 0.5, - 0.5, - r/2);
	girdle[7] = new Point3D(- lw * 0.5,   0.5, - r/2);
}
