struct node
{
	int num;
	struct node *next;
};

struct ptrList
{
	struct node* currentPtr;
	struct ptrList* nextList;
};

struct coords
{
	double x;
	double y;
	double z;
	double nf;
	struct coords* next;
};

extern void add_to_node(struct node* root, int el, int*numAdded);
extern void free_node(struct node* ptr);
extern void free_ptrList(struct ptrList *ptr, int NumEl);
extern int add_to_coords(struct coords* root, double x_add, double y_add, double z_add, double nf_add, int* numAdded);
