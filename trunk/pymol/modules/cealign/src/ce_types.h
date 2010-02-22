/*
// Typical XYZ point and array of points
*/
typedef struct {
	double x;
	double y;
	double z;
} cePoint, *pcePoint;

/*
// An AFP (aligned fragment pair), and list/pointer
*/
typedef struct {
	int first;
	int second;
} afp, *path, **pathCache;
