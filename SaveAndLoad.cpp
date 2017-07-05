template<typename Unit>
int get_position(list<Unit> *unit_list, typename list<Unit>::iterator unit_iter)
{
	//cout << "i'm here" << endl;
	//cout << "size:" << unit_list->size() << endl;
	//system("pause");
	int position = 0;
	if (unit_list->size() != 0 && unit_iter != unit_list->end())
		//if (false)
	{
		//cout << "lalala" << endl;
		//system("pause");
		auto iter = unit_list->begin();

		while (iter != unit_iter)
		{
			//cout << "position:" << position << endl;
			//system("pause");
			/*if (position>1000)
			{
				cout << "comparison error" << endl;
				system("pause");
			}*/
			position++;
			iter++;
		}
	}
	else
	{
		//cout << "lululu" << endl;
		//system("pause");
		position = -1;
	}
	//cout << "position:" << position << endl;
	return position;
};
