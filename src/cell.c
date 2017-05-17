#include "cell.h"

#include "space.h"
#include "potential.h"
#include "mympi.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// 初始化细胞链表
void initCells(struct SpacialStr* space, struct PotentialStr* potential, struct CellStr** cel){

	*cel = (Cell*)malloc(sizeof(Cell));
  Cell* cells = *cel;

	// 保证细胞长度大于等于截断距离
	for (int i = 0; i < 3; i++)
   	{
      	cells->xyzCellNum[i] = space->myLength[i] / potential->cutoff; 
      	cells->cellLength[i] = space->myLength[i] / ((double) cells->xyzCellNum[i]);
   	}

   	// 实际细胞数为 x * y * z
   	cells->myCellNum = cells->xyzCellNum[0] * cells->xyzCellNum[1] * cells->xyzCellNum[2];
   
   	// 通信细胞数
   	cells->commCellNum = 2 * ((cells->xyzCellNum[0] + 2) *
                         (cells->xyzCellNum[1] + cells->xyzCellNum[2] + 2) +
                         (cells->xyzCellNum[1] * cells->xyzCellNum[2]));

   	// 总细胞数
   	cells->totalCellNum = cells->myCellNum + cells->commCellNum;
   
   	cells->atomNum = malloc(cells->totalCellNum*sizeof(int));

   	// 初始化各细胞中原子数为0
   	for (int i = 0; i < cells->totalCellNum; i++)
      	cells->atomNum[i] = 0;

    ///test 测试cell坐标正确性

    // if(getMyRank()== 5){
    //     int3 xyz;
    //     int j;
    //     for(int i = 0; i < cells->totalCellNum; i++)
    //     {
    //         getXYZByCell(cells,xyz, i);
    //         printf("i: %d xyz:%d %d %d\n",i,xyz[0],xyz[1],xyz[2]);
    //         j = findCellByXYZ(cells, xyz);
    //         printf("j: %d xyz:%d %d %d\n\n",j,xyz[0],xyz[1],xyz[2]);
    //     }
    // }
}

// 根据原子坐标找到所在的细胞
int findCellByCoord(Cell* cells, Spacial* space, double3 coord){

    double* myMin = space->myMin;
    double* myMax = space->myMax;
    int*    xyzCellNum = cells->xyzCellNum; 

    // 细胞所在的位置
    int3 cellPos; 

    cellPos[0] = (int)(floor((coord[0] - myMin[0])/cells->cellLength[0]));
    cellPos[1] = (int)(floor((coord[1] - myMin[1])/cells->cellLength[1]));
    cellPos[2] = (int)(floor((coord[2] - myMin[2])/cells->cellLength[2]));

    // 如果原子坐标超出了空间边界，则加入至通信区域的细胞中
    for(int i = 0; i< 3 ; i++){
        if(coord[i] >= myMax[i])
            cellPos[i] = xyzCellNum[i];
    }

    return findCellByXYZ(cells, cellPos);
}

// 根据细胞位置xyz返回细胞序号，即该空间中第几个细胞
int findCellByXYZ(Cell* cells, int* xyz){

    int cell;

    int myCellNum = cells->myCellNum;
    int *xyzCellNum = cells->xyzCellNum;

    // Z轴正方向的通信区域细胞
    if (xyz[2] == xyzCellNum[2])
        cell = myCellNum + 2*xyzCellNum[2]*xyzCellNum[1] + 2*xyzCellNum[2]*(xyzCellNum[0]+2) +
            (xyzCellNum[0]+2)*(xyzCellNum[1]+2) + (xyzCellNum[0]+2)*(xyz[1]+1) + (xyz[0]+1);
    // Z轴负方向的通信区域细胞
    else if (xyz[2] == -1)
        cell = myCellNum + 2*xyzCellNum[2]*xyzCellNum[1] + 2*xyzCellNum[2]*(xyzCellNum[0]+2) +
            (xyzCellNum[0]+2)*(xyz[1]+1) + (xyz[0]+1);
    // Y轴正方向的通信区域细胞
    else if (xyz[1] == xyzCellNum[1])
        cell = myCellNum + 2*xyzCellNum[2]*xyzCellNum[1] + xyzCellNum[2]*(xyzCellNum[0]+2) +
            (xyzCellNum[0]+2)*xyz[2] + (xyz[0]+1);
   // Y轴负方向的通信区域细胞
    else if (xyz[1] == -1)
        cell = myCellNum + 2*xyzCellNum[2]*xyzCellNum[1] + xyz[2]*(xyzCellNum[0]+2) + (xyz[0]+1);
   // X轴正方向的通信区域细胞
    else if (xyz[0] == xyzCellNum[0])
        cell = myCellNum + xyzCellNum[1]*xyzCellNum[2] + xyz[2]*xyzCellNum[1] + xyz[1];
   // X轴负方向的通信区域细胞
    else if (xyz[0] == -1)
        cell = myCellNum + xyz[2]*xyzCellNum[1] + xyz[1];
   // 本空间中实际的细胞
    else
        cell = xyz[0] + xyzCellNum[0]*xyz[1] + xyzCellNum[0]*xyzCellNum[1]*xyz[2];

    return cell;
}

// 根据细胞序号返回细胞位置xyz,与函数findCellByXYZ互为逆过程
void getXYZByCell(Cell* cells,int *xyz, int num){

    int *xyzCellNum = cells->xyzCellNum;
   
    if( num < cells->myCellNum)
    {
        xyz[0] = num % xyzCellNum[0];
        num /= xyzCellNum[0];
        xyz[1] = num % xyzCellNum[1];
        xyz[2] = num / xyzCellNum[1];
    }
    else 
    {
        int ink;
        ink = num - cells->myCellNum;
        if (ink < 2*xyzCellNum[1]*xyzCellNum[2])
        {
            if (ink < xyzCellNum[1]*xyzCellNum[2]) 
            {
                xyz[0] = 0;
            }
            else 
            {
                ink -= xyzCellNum[1]*xyzCellNum[2];
                xyz[0] = xyzCellNum[0] + 1;
            }
            xyz[1] = 1 + ink % xyzCellNum[1];
            xyz[2] = 1 + ink / xyzCellNum[1];
        }
        else if (ink < (2 * xyzCellNum[2] * (xyzCellNum[1] + xyzCellNum[0] + 2))) 
        {
            ink -= 2 * xyzCellNum[2] * xyzCellNum[1];
            if (ink < ((xyzCellNum[0] + 2) *xyzCellNum[2])) 
            {
                xyz[1] = 0;
            }
            else 
            {
                ink -= (xyzCellNum[0] + 2) * xyzCellNum[2];
                xyz[1] = xyzCellNum[1] + 1;
            }
            xyz[0] = ink % (xyzCellNum[0] + 2);
            xyz[2] = 1 + ink / (xyzCellNum[0] + 2);
        }
        else 
        {
            ink -= 2 * xyzCellNum[2] * (xyzCellNum[1] + xyzCellNum[0] + 2);
            if (ink < ((xyzCellNum[0] + 2) * (xyzCellNum[1] + 2))) 
            {
                xyz[2] = 0;
            }
            else 
            {
                ink -= (xyzCellNum[0] + 2) * (xyzCellNum[1] + 2);
                xyz[2] = xyzCellNum[2] + 1;
            }
            xyz[0] = ink % (xyzCellNum[0] + 2);
            xyz[1] = ink / (xyzCellNum[0] + 2);
        }
        xyz[0]--;
        xyz[1]--;
        xyz[2]--;
    }
}

// // 根据细胞位置xyz,返回在共享内存中的细胞序号,若不是共享内存内,则返回-1
// int getSMCellByXYZ(Cell* cells, int* xyz){

//     int *xyzCellNum = cells->xyzCellNum;
//     //如果不在共享内存范围内
//     if(xyz[0]>0 && (xyz[0]<xyzCellNum[0]-1)&&xyz[1]>0 && (xyz[1]<xyzCellNum[1]-1)&&xyz[2]>0 && (xyz[2]<xyzCellNum[2]-1)){
//         return -1;   
//     }
// }     