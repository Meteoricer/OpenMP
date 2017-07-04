template<typename Unit>
int get_position(list<Unit> unit_list, typename list<Unit>::iterator unit_iter)
{
	int position = 0;
	if (unit_list.size() != 0 && unit_iter != unit_list.end())
	{
		auto iter = unit_list.begin();

		while (iter != unit_iter)
		{
			position++;
			iter++;
		}
	}
	else
	{
		position = -1;
	}
	return position;
};
