#include <stdlib.h>
#include <math.h>
#include "Data_Structures.h"

void add_to_node(struct node* root, int el, int* numAdded)
{
	struct node* ptr = root;
	int cont = 1;
	int add = 1;
	while (cont)
	{
		int num = ptr->num;
		if (el == num)
		{
			add = 0;
			cont = 0;
		}
		if (num == -1)
		{
			if (add)
			{
				ptr->num = el;
				ptr->next = (struct node*) malloc(sizeof(struct node));
				ptr = ptr->next;
				ptr->num = -1;

				(*numAdded)++;

			}
			cont = 0;
		}
		else
			ptr = ptr->next;
	}

}

int add_to_coords(struct coords* root, double x_add, double y_add, double z_add, double nf_add, int* numAdded)
{
	struct coords* ptr = root;
	int cont = 1;
	int add = 1;
	int nodeNum = -1;
	int count = 0;
	if (ptr == NULL)
	{
		add = 1;
		cont = 0;
	}
	else
	{	
		while (ptr != NULL && cont)
		{
			double m_x = ptr->x;
			double m_y = ptr->y;
			if ((fabs(m_x - x_add) < 1e-8) && (fabs (m_y - y_add) < 1e-8))
			{
				add = 0;
				cont = 0;
				nodeNum = count;
			}
			ptr = ptr->next;
			count ++;
		}
	}

	if (add)
	{
		ptr = (struct coords*) malloc(sizeof(struct coords));
		ptr->x = x_add;
		ptr->y = y_add;
		ptr->z = z_add;
		ptr->nf = nf_add;
		ptr->next = NULL;
		nodeNum = *numAdded;
		(*numAdded)++;
	}

	return nodeNum;
}

void free_node(struct node *ptr)
{
	while(ptr->num != -1)
	{
		struct node *currPtr = ptr;
		ptr = ptr->next;
		free(currPtr);
	}
	free(ptr);
}

void free_ptrList(struct ptrList *ptr, int NumEl)
{
	for (int i = 0; i < NumEl; i++)
	{
		struct ptrList *currPtrList = ptr;
		ptr = ptr->nextList;
		free_node(currPtrList->currentPtr);
		free(currPtrList);
	}
}
